#ifndef TEMPLAT_FFT_EXTERNAL_PARAFAFT_PARAFAFTPLANHOLDER_H
#define TEMPLAT_FFT_EXTERNAL_PARAFAFT_PARAFAFTPLANHOLDER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2026

#ifndef NOFFT
#ifdef HAVE_MPI
#ifdef HAVE_PARAFAFT
#include <parafaft_r2c.hpp>
#endif
#endif
#endif

#include "TempLat/util/exception.h"
#include "TempLat/fft/fftlibraryinterface.h"
#include "TempLat/parallel/mpi/cartesian/mpicartesiangroup.h"
#include "TempLat/lattice/memory/memoryblock.h"

#include <memory>
#include "TempLat/parallel/device.h"

namespace TempLat
{
  // Depending on the device, we alias ParaFaFT_Backend to the appropriate backend type. Parafaft supports multiple
  // backends (FFTW, cuFFT, HIPFFT), and we select the one based on compile-time definitions set by CMake when detecting
  // the device and available libraries. This allows us to write device-agnostic code in the ParafaftPlanHolder, while
  // still leveraging the performance benefits of the appropriate backend for the target platform. Each backend is
  // itself a class template parameterized on the floating-point precision, so the alias is a template alias.
#ifdef DEVICE_CUDA
  template <typename T> using ParaFaFT_Backend = parafaft::CuFFTBackend<T>;
#elif defined(CL_HIP)
  template <typename T> using ParaFaFT_Backend = parafaft::HipFFTBackend<T>;
#else
  template <typename T> using ParaFaFT_Backend = parafaft::FFTWBackend<T>;
#endif

  MakeException(ParafaftPlanHolderException);
  MakeException(ParafaftCompiledWithoutSinglePrecisionSupport);

  /** @brief Plan holder for parafaft FFT transforms.
   *
   * This class holds a parafaft::ParaFaFT_R2C object and implements the
   * PlanInterface<T> methods for r2c and c2r transforms.
   *
   * Unit test: ctest -R test-parafaftplanholder
   **/
  template <typename T, size_t NDim> class ParafaftPlanHolder : public FFTPlanInterface<T, NDim>
  {
  public:
    using Complex = complex<T>;

#ifdef HAVE_MPI
#ifdef HAVE_PARAFAFT
    /**
     * @brief Constructor - takes ownership of the parafaft object.
     *
     * We store the MPICartesianGroup to keep the MPI communicator alive.
     */
    ParafaftPlanHolder(MPICartesianGroup group,
                       std::shared_ptr<parafaft::ParaFaFT_R2C<NDim, ParaFaFT_Backend<T>>> parafaftObj)
        : mGroup(group), mParafaft(parafaftObj)
    {
    }
#endif
#endif

    virtual ~ParafaftPlanHolder()
    {
      // Parafaft cleans up internally via its destructor
    }

    /**
     * @brief Complex-to-real backward transform (in-place).
     *
     * Parafaft's backward() copies input to internal buffers before processing,
     * so using the same memory for both is safe.
     */
    virtual void c2r(MemoryBlock<T, NDim> &mBlock) override { execute_c2r(mBlock); }

    /**
     * @brief Real-to-complex forward transform (in-place).
     *
     * Parafaft's forward() copies input to internal buffers before processing,
     * so using the same memory for both is safe.
     */
    virtual void r2c(MemoryBlock<T, NDim> &mBlock) override { execute_r2c(mBlock); }

  private:
#ifdef HAVE_MPI
#ifdef HAVE_PARAFAFT
    // Keep group alive for MPI communicator lifetime
    MPICartesianGroup mGroup;

    // Shared pointer to parafaft object
    std::shared_ptr<parafaft::ParaFaFT_R2C<NDim, ParaFaFT_Backend<T>>> mParafaft;

#endif
#endif

    void execute_r2c(MemoryBlock<T, NDim> &mBlock)
    {
#ifdef HAVE_MPI
#ifdef HAVE_PARAFAFT
      // Parafaft's forward_in_place accepts padded buffers directly
      // Buffer layout: [N0_local][N1_local][2*(N/2+1)] - matches CosmoLattice
      mParafaft->forward_in_place(mBlock.data());
#endif
#endif
    }

    void execute_c2r(MemoryBlock<T, NDim> &mBlock)
    {
#ifdef HAVE_MPI
#ifdef HAVE_PARAFAFT
      // Parafaft's backward_in_place accepts padded buffers directly
      // Buffer layout: [N0_local][N1_local][2*(N/2+1)] - matches CosmoLattice
      mParafaft->backward_in_place(mBlock.data());
#endif
#endif
    }
  };
} // namespace TempLat

#endif
