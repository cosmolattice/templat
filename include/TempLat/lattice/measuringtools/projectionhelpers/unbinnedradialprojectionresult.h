#ifndef TEMPLAT_LATTICE_MEASUREMENTS_PROJECTIONHELPERS_UNBINNEDRADIALPROJECTIONRESULT_H
#define TEMPLAT_LATTICE_MEASUREMENTS_PROJECTIONHELPERS_UNBINNEDRADIALPROJECTIONRESULT_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros, Year: 2026
//            Based on: Wessel Valkenburg, Year: 2019

#include "TempLat/util/exception.h"
#include "TempLat/lattice/algebra/helpers/getfloattype.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionsinglebinandvalue.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionsinglequantity.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionrebinner.h"

#include "TempLat/parallel/device.h"

namespace TempLat
{

  // Own exception so this header is self-contained (radialprojectionresult.h declares a separate one).
  MakeException(UnbinnedRadialProjectionResultFinalizationException);

  /** @brief A class which holds the result of a radial projection without binning.
   *  For each possible value of the Fourier coordinate it saves:
   *  the average value,
   *  the sample variance (= <f^2> - <f>^2 ) of the values,
   *  and the number of entries in a bin (multiplicity).
   *
   *  Upon call to finalize, constructs the parent vector, which you can use as you wish.
   *  This contains pairs for each value of the momentum, with the momentum (first) and
   *  the value result (second), as a `RadialProjectionSingleDatum'
   *
   *  Unit test: ctest -R test-unbinnedradialprojectionresult
   **/

  template <typename T>
  class UnbinnedRadialProjectionResult : public std::vector<std::tuple<T, RadialProjectionSingleDatum<T>>>
  {
  public:
    using floatType = typename GetFloatType<T>::type;

    // Put public methods here. These should change very little over time.
    UnbinnedRadialProjectionResult(size_t nBins, bool pIsInFourier = false)
        : std::vector<std::tuple<T, RadialProjectionSingleDatum<T>>>(), finalizedOnce(false), mNBins(nBins),
          mValues(mNBins), mIsInFourier(pIsInFourier)
    {
      mMultiplicitiesDevice = DeviceView("RadialProjectionResult::mMultiplicitiesDevice", mNBins);
      mMultiplicities = device::memory::createMirrorView(mMultiplicitiesDevice);
    }

    /** \brief Rescale the results with a function of x or k (bin location),
     *  using for now a simple lambda function of single float parameter, which
     *  is your x in f(x). In other words, you give f(x).
     */
    template <typename LL> UnbinnedRadialProjectionResult &rescale(LL rescaler)
    {
      // for (auto&& it : *this) {
      for (size_t i = 0; i < this->size(); ++i) {
        std::get<1>((*this)[i]) *= rescaler(std::get<0>((*this)[i]));
      }
      return *this;
    }

    /** @brief Rescale the bin positions with a normalization (for example dimensionful).
     */
    UnbinnedRadialProjectionResult &rescaleBins(T scale)
    {
      for (auto &&it : *this) {
        std::get<0>(it) *= scale;
      }
      return *this;
    }

    auto integrate(bool dummy)
    {
      auto total = 0.;
      for (size_t i = 0; i < this->size(); ++i) {
        total += std::get<1>((*this)[i]).average / std::get<0>((*this)[i])  * binWidths[i];
      }
      return total;
    }

    UnbinnedRadialProjectionResult &sumInsteadOfAverage()
    {
      floatType intMultiplicity = 0;

      for (auto &&it : *this) {
        intMultiplicity =
            std::get<1>(it).multiplicity *
            (mIsInFourier ? 2 : 1); // Multiply by 2 if in Fourier space (only half of the last coordinate is iterated
        // over in this case because of reflection symmetry for real data).
        std::get<1>(it).average *= intMultiplicity;
      }
      return *this;
    }

    UnbinnedRadialProjectionResult &renormalizeBins()
    {
      /*After finalizing, we scale the bins so that they can be compared to the binned power psectrum. This is equal to
       * considering the unbinned power spectrum with bins of variable width*/
      setBinWidths();

      const size_t n = (*this).size();
      if (n < 2) return *this; // need at least two bins to form a finite bin width

      // Divide a bin's average by its width, guarding against coincident positions (width 0).
      auto rescale = [&](size_t i, floatType width) {
        if (width != floatType(0)) {
          std::get<1>((*this)[i]).average /= width;
        }
      };

      for (size_t i = 0; i < n; ++i) rescale(i, binWidths[i]);

      return *this;
    }

    void setBinWidths()
    {
      /*After finalizing, we scale the bins so that they can be compared to the binned power psectrum. This is equal to
       * considering the unbinned power spectrum with bins of variable width*/
      const size_t n = (*this).size();
      if (n < 2) return; // need at least two bins to form a finite bin width


      // First bin: forward half-difference, minus the 0.5 offset of the lowest binned mode.
      binWidths[0] = (std::get<0>((*this)[1]) - std::get<0>((*this)[0])) / 2 - 0.5;
      // Interior bins: central difference.
      for (size_t i = 1; i + 1 < n; ++i) {
        binWidths[i] = (std::get<0>((*this)[i + 1]) - std::get<0>((*this)[i - 1])) / 2;
      }
      // Last bin: backward difference (was an out-of-bounds read/write at index size()).
      binWidths[n-1] = std::get<0>((*this)[n - 1]) - std::get<0>((*this)[n - 2]);
    }

    std::string toString(int verbosity = 0) const
    {
      if ((device::Idx)this->size() < 1) return "";
      std::stringstream sstream;
      // Each entry is a (bin position, datum) tuple; the datum carries the header/value formatting.
      sstream << std::get<1>(this->front()).getHeader(1, "#", true, verbosity) << "\n";
      for (auto &&it : *this) {
        sstream << std::get<0>(it) << " " << std::get<1>(it).toString(true, verbosity) << "\n";
      }
      return sstream.str();
    }

    friend std::ostream &operator<<(std::ostream &ostream, const UnbinnedRadialProjectionResult &rpr)
    {
      ostream << rpr.toString() << "\n";
      return ostream;
    }

    template <typename S> friend class RadialProjector;

    template <typename S>
    friend UnbinnedRadialProjectionResult<S> operator+(const UnbinnedRadialProjectionResult<S> &a,
                                                       const UnbinnedRadialProjectionResult<S> &b);

    DEVICE_FORCEINLINE_FUNCTION
    void add_device(device::Idx i, const T &value, const T &weight = (T)1) const
    {
      mValues.add_device(i, value, weight);
      device::atomic_add(&mMultiplicitiesDevice(i), weight);
    }

    /** @brief RadialProjector calls this as the last step, does the transposition of the result vecotrs into one vector
     * of RadialProjectionSingleBinAndValue<T>. */
    UnbinnedRadialProjectionResult &finalize(MPICommReference comm)
    {
      if (finalizedOnce)
        throw UnbinnedRadialProjectionResultFinalizationException("Can only finalize once per instance.");

      finalizedOnce = true;
      pull();

      comm.Allreduce(mMultiplicities, MPI_SUM);

      mValues.finalize(comm);

      auto &mFullResult = *this;

      mFullResult.clear();
      for (size_t i = 0; i < mNBins; ++i) {
        if (mMultiplicities[i] > 0.) {
          std::tuple<T, RadialProjectionSingleDatum<T>> next(static_cast<T>(std::sqrt(i)), mValues.getFinal(i, mMultiplicities[i]));
          this->push_back(next);
        }
      }
      binWidths.resize(this->size(), 1.);

      return *this;
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    bool finalizedOnce;
    size_t mNBins;
    RadialProjectionSingleQuantity<T> mValues;
    std::vector<T> binWidths; // Width of the bins

    using DeviceView = device::memory::NDView<floatType, 1>;
    using HostMirror = typename DeviceView::host_mirror_type;

    HostMirror mMultiplicities;
    DeviceView mMultiplicitiesDevice;

    bool mIsInFourier;

    void pull()
    {
      mValues.pull();
      device::memory::copyDeviceToHost(mMultiplicitiesDevice, mMultiplicities.data());
    }

    void push()
    {
      mValues.push();
      device::memory::copyHostToDevice(mMultiplicities.data(), mMultiplicitiesDevice);
    }
  };

  template <typename T>
  UnbinnedRadialProjectionResult<T> operator+(const UnbinnedRadialProjectionResult<T> &a,
                                              const UnbinnedRadialProjectionResult<T> &b)
  {
    UnbinnedRadialProjectionResult<T> res(a);

    for (size_t i = 0; i < a.size(); ++i) {
      std::get<1>(res[i]).average += std::get<1>(b[i]).average;
    }
    return res;
  }

  template <typename T, class R>
  UnbinnedRadialProjectionResult<T> operator*(R &&func, const UnbinnedRadialProjectionResult<T> &obj)
  {
    UnbinnedRadialProjectionResult<T> res(obj);
    return res.rescale(func);
  }

  template <typename T>
  UnbinnedRadialProjectionResult<T> operator*(double scale, const UnbinnedRadialProjectionResult<T> &obj)
  {
    auto func = [&](auto x) { return scale; };
    return func * obj;
  }

  template <typename T>
  UnbinnedRadialProjectionResult<T> operator*(float scale, const UnbinnedRadialProjectionResult<T> &obj)
  {
    auto func = [&](auto x) { return scale; };
    return func * obj;
  }

  template <typename T>
  UnbinnedRadialProjectionResult<T> operator*(int scale, const UnbinnedRadialProjectionResult<T> &obj)
  {
    auto func = [&](auto x) { return scale; };
    return func * obj;
  }
} // namespace TempLat

#endif
