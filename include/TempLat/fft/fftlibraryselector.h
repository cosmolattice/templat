#ifndef TEMPLAT_FFT_FFTLIBRARYSELECTOR_H
#define TEMPLAT_FFT_FFTLIBRARYSELECTOR_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Franz R. Sattler,  Year: 2025

#include "TempLat/fft/external/fftw/fftwinterface.h"
#include "TempLat/fft/fftdecomposition.h"

#ifndef NOFFT

#ifdef HAVE_MPI
#ifdef HAVE_PARAFAFT
#include "TempLat/fft/external/parafaft/parafaftinterface.h"
#endif
#endif

#ifdef HAVE_KOKKOSFFT
#include "TempLat/fft/external/kokkosfft/kokkosfftinterface.h"
#endif

#endif

namespace TempLat
{
  MakeException(FFTLibraryDoubleInitializationException);
  MakeException(FFTLibraryDecompositionMismatchException);

  /**
   * @brief I wrapped this in a struct for a very specific case: If we have multiple translation units (cpp files)
   * which include this header, and each calls getFFTSessionGuards, then we will have multiple static variables, one
   * per translation unit, and the guards will not work as intended. By wrapping it in a struct, we ensure that
   * there is only one instance of the static variable, no matter how many translation units include this header.
   *
   */
  struct holdStaticGuard {
    /** @brief An inline function for storing a static global variable in a header. A lock that
        verifies that we do not accidentally call the FFT initi/fin-alizations twice. */
    static bool getSessionGuardsWasCalledOnce()
    {
      static bool wasOnce = false;
      bool result = wasOnce;
      wasOnce = true;
      return result;
    }
  };

  static inline std::vector<std::shared_ptr<FFTSessionGuard>> getFFTSessionGuards(bool pVerbose = true)
  {
    if (holdStaticGuard::getSessionGuardsWasCalledOnce())
      throw FFTLibraryDoubleInitializationException("You can only call getSessionGuards once.");

    std::vector<std::shared_ptr<FFTSessionGuard>> result;

    /* add your guards here. */
    /* Note: the standard guarantees that the destructors are called in the
     inverse order of appearance in the vector. So if you must construct after FFTW and
     destruct before FFTW, you should be safe if you add your thing after FFTW.
     */
#ifndef NOFFT
    result.push_back(getFFTWSessionGuard(pVerbose));
#ifdef HAVE_MPI

#endif
#endif

    return result;
  }

  /** @brief A class which sets up the interface with the appropriate FFT library.
   * Once you have implemented the FFTLibraryInterface for your library, add it to the logic here.
   *
   * Unit test: ctest -R test-fftlibraryselector
   **/
  template <size_t NDim> class FFTLibrarySelector
  {
  public:
    // Put public methods here. These should change very little over time.
    FFTLibrarySelector(MPICartesianGroup group, const device::IdxArray<NDim> &nGridPoints,
                       bool forbidTransposition = false)
        : mGroup(group), mNGridPoints(nGridPoints), theLibrary(nullptr), mLayout(mNGridPoints), madePlansFloat(false),
          madePlansDouble(false), verbose(false)
    {
      std::tie(theLibrary, backend) = selectBackend(group.size());
      if (mGroup.getRank() == 0) sayShort << "Using " << backend << " backend for FFTs.\n";
      verifyDecompositionMatchesBackend(group, nGridPoints);
      mLayout = theLibrary->computeLocalSizes(mGroup, mNGridPoints, forbidTransposition);
    }

    const auto &getLayout() { return mLayout; }

    const std::string &getBackend() const { return backend; }

    void setVerbose() { verbose = true; }

    /** @brief The decomposition the runtime-selected FFT backend will use for this setup.
     *
     *  Used by `FFTMPIDomainSplit` to build the MPI Cartesian group so that the group and the
     *  backend always agree. The priority chain here must stay in lockstep with `selectBackend`
     *  below — both call the same set of backends in the same order under the same `#ifdef`s,
     *  so the pre-topology query and the post-topology selection cannot pick different backends.
     */
    static FFTDecomposition<NDim> decomposition(MPICommReference baseComm,
                                                const device::IdxArray<NDim> &nGridPoints)
    {
      const ptrdiff_t nProcesses = baseComm.size();
#ifdef HAVE_MPI
#ifdef HAVE_PARAFAFT
      if constexpr (NDim >= 2) {
        if (nProcesses > 1) return ParafaftInterface<NDim>::decomposition(baseComm, nGridPoints);
      }
#endif
#endif
#ifdef HAVE_KOKKOSFFT
      if (nProcesses == 1) {
        if constexpr (NDim <= 3) return KokkosFFTInterface<NDim>::decomposition(baseComm, nGridPoints);
      }
#endif
      (void) nProcesses;
      return FFTWInterface<NDim>::decomposition(baseComm, nGridPoints);
    }

    void r2c(MemoryBlock<double, NDim> &mBlock)
    {
      getPlans_double();
      if (verbose) sayMPI << "Going to perform double r2c.\n";
      mPlansDouble->r2c(mBlock);
      mBlock.flagHostMirrorOutdated();
    }

    void r2c(MemoryBlock<float, NDim> &mBlock)
    {
      getPlans_float();
      if (verbose) sayMPI << "Going to perform float r2c.\n";
      mPlansFloat->r2c(mBlock);
      mBlock.flagHostMirrorOutdated();
    }

    void c2r(MemoryBlock<double, NDim> &mBlock)
    {
      getPlans_double();
      if (verbose) sayMPI << "Going to perform double c2r.\n";
      mPlansDouble->c2r(mBlock);
      mBlock.flagHostMirrorOutdated();
    }

    void c2r(MemoryBlock<float, NDim> &mBlock)
    {
      getPlans_float();
      if (verbose) sayMPI << "Going to perform float c2r.\n";
      mPlansFloat->c2r(mBlock);
      mBlock.flagHostMirrorOutdated();
    }

    void getPlans_float()
    {
      if (!madePlansFloat) {
        if (verbose) sayMPI << "Going to prepare float FFT plans.\n";
        madePlansFloat = true;
        mPlansFloat = theLibrary->getPlans_float(mGroup, mLayout);
      }
    }
    void getPlans_double()
    {
      if (!madePlansDouble) {
        if (verbose) sayMPI << "Going to prepare double FFT plans.\n";
        madePlansDouble = true;
        mPlansDouble = theLibrary->getPlans_double(mGroup, mLayout);
      }
    }

  private:
    /** @brief Reject a user-supplied `MPICartesianGroup` whose shape the selected backend cannot use.
     *
     *  Backends fall into two classes:
     *  - No explicit grid (`dims` all zero): only `nDimsToSplit` is constrained — the rank count
     *    must be factored over at most that many leading dimensions. Groups that split more
     *    dimensions (or, for FFTW, any secondary dimension) are rejected.
     *  - Explicit grid: every entry of the provided decomposition must match the probed dims.
     *  Auto-built groups from `FFTMPIDomainSplit::makeMPIGroup` always pass.
     */
    void verifyDecompositionMatchesBackend(const MPICartesianGroup &group,
                                           const device::IdxArray<NDim> &nGridPoints) const
    {
      const FFTDecomposition<NDim> expected = decomposition(group.getBaseComm(), nGridPoints);
      const std::vector<int> &actual = group.getDecomposition();
      if (actual.size() != NDim)
        throw FFTLibraryDecompositionMismatchException("Cartesian group decomposition has ", actual.size(),
                                                       " entries, expected ", NDim, ".");

      bool explicitGrid = expected.nDimsToSplit > 0;
      for (ptrdiff_t i = 0; i < expected.nDimsToSplit; ++i)
        if (expected.dims[i] <= 0) { explicitGrid = false; break; }

      if (explicitGrid) {
        for (size_t i = 0; i < NDim; ++i) {
          if (actual[i] != expected.dims[i])
            throw FFTLibraryDecompositionMismatchException(
                "FFT backend ", backend, " requires MPI decomposition {",
                expectedDimsToString(expected), "} but the provided MPICartesianGroup has {",
                actualDimsToString(actual),
                "}. Use FFTMPIDomainSplit::makeMPIGroup to build a matching group.");
        }
      } else {
        ptrdiff_t actualSplits = 0;
        for (size_t i = 0; i < NDim; ++i)
          if (actual[i] > 1) ++actualSplits;
        if (actualSplits > expected.nDimsToSplit)
          throw FFTLibraryDecompositionMismatchException(
              "FFT backend ", backend, " supports splitting at most ", expected.nDimsToSplit,
              " dimension(s) but the provided MPICartesianGroup splits ", actualSplits,
              ". Use FFTMPIDomainSplit::makeMPIGroup to build a matching group.");
      }
    }

    static std::string expectedDimsToString(const FFTDecomposition<NDim> &d)
    {
      std::string s;
      for (size_t i = 0; i < NDim; ++i) {
        if (i) s += ",";
        s += std::to_string(d.dims[i]);
      }
      return s;
    }
    static std::string actualDimsToString(const std::vector<int> &v)
    {
      std::string s;
      for (size_t i = 0; i < v.size(); ++i) {
        if (i) s += ",";
        s += std::to_string(v[i]);
      }
      return s;
    }

    /** @brief Shared backend-selection priority chain. Structure must stay in lockstep with the
     *  `decomposition` static above so that the decomposition query and the runtime backend pick
     *  never disagree. Any change to the order or predicates here must be mirrored there.
     */
    static std::pair<std::shared_ptr<FFTLibraryInterface<NDim>>, std::string> selectBackend(ptrdiff_t nProcesses)
    {
#ifdef HAVE_MPI
#ifdef HAVE_PARAFAFT
      if constexpr (NDim >= 2) {
        if (nProcesses > 1) return {std::make_shared<ParafaftInterface<NDim>>(), "Parafaft"};
      }
#endif
#endif
#ifdef HAVE_KOKKOSFFT
      if (nProcesses == 1) {
        if constexpr (NDim <= 3) return {std::make_shared<KokkosFFTInterface<NDim>>(), "KokkosFFT"};
      }
#endif
      (void) nProcesses;
      return {std::make_shared<FFTWInterface<NDim>>(), "FFTW"};
    }

    /* Put all member variables and private methods here. These may change arbitrarily. */
    MPICartesianGroup mGroup;
    device::IdxArray<NDim> mNGridPoints;
    std::shared_ptr<FFTLibraryInterface<NDim>> theLibrary;
    FFTLayoutStruct<NDim> mLayout;

    bool madePlansFloat;
    bool madePlansDouble;
    std::shared_ptr<FFTPlanInterface<float, NDim>> mPlansFloat;
    std::shared_ptr<FFTPlanInterface<double, NDim>> mPlansDouble;

    bool verbose;

    std::string backend;
  };
} // namespace TempLat

#endif
