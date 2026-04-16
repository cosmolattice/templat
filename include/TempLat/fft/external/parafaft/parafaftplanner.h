#ifndef TEMPLAT_FFT_EXTERNAL_PARAFAFT_PARAFAFTPLANNER_H
#define TEMPLAT_FFT_EXTERNAL_PARAFAFT_PARAFAFTPLANNER_H

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
#include "TempLat/fft/external/parafaft/parafaftplanholder.h"
#include "TempLat/lattice/memory/memoryblock.h"

#include <memory>

namespace TempLat
{
  MakeException(ParafaftPlannerException);

  /** @brief Plan creation for parafaft FFT transforms.
   *
   * Creates parafaft::ParaFaFT_R2C objects wrapped in ParafaftPlanHolder.
   * Supports both double and (when HAVE_FFTFLOAT is defined) single precision.
   *
   * Unit test: ctest -R test-parafaftplanner
   **/
  template <size_t NDim> class ParafaftPlanner : public FFTLibraryInterface<NDim>
  {
  public:
    ParafaftPlanner() {}

    /**
     * @brief Set planner patience level.
     *
     * Parafaft uses FFTW internally with FFTW_ESTIMATE, so this is a no-op.
     * The parameter is accepted for interface compatibility.
     */
    virtual void setPlannerPatience(int level) override
    {
      // Parafaft uses FFTW_ESTIMATE internally - no patience setting
      (void)level;
    }

    /**
     * @brief Create float-precision plans.
     *
     * Enabled when HAVE_FFTFLOAT is defined (CMake FLOAT=ON). Requires libfftw3f
     * on CPU/FFTW builds (parafaft's PARAFAFT_FFTW3F_AVAILABLE); TempLat's CMake
     * already enforces this at configure time when FLOAT=ON.
     */
    virtual std::shared_ptr<FFTPlanInterface<float, NDim>> getPlans_float(const MPICartesianGroup &group,
                                                                          const FFTLayoutStruct<NDim> &layout) override
    {
#ifndef HAVE_FFTFLOAT
      throw ParafaftCompiledWithoutSinglePrecisionSupport("CosmoLattice compiled without float FFT support.");
      return std::shared_ptr<FFTPlanInterface<float, NDim>>();
#elif !defined(HAVE_MPI) || !defined(HAVE_PARAFAFT)
      throw ParafaftPlannerException("Parafaft requires MPI and HAVE_PARAFAFT.");
      return std::shared_ptr<FFTPlanInterface<float, NDim>>();
#else
      return make_plans<float>(group, layout);
#endif
    }

    /**
     * @brief Create double-precision plans.
     *
     * Creates a parafaft::ParaFaFT_R2C<NDim, ParaFaFT_Backend<double>> object and wraps it in ParafaftPlanHolder.
     * Parafaft's in-place API accepts padded buffers matching CosmoLattice's layout.
     */
    virtual std::shared_ptr<FFTPlanInterface<double, NDim>>
    getPlans_double(const MPICartesianGroup &group, const FFTLayoutStruct<NDim> &layout) override
    {
#if defined(HAVE_MPI) && defined(HAVE_PARAFAFT)
      return make_plans<double>(group, layout);
#else
      throw ParafaftPlannerException("Parafaft requires MPI and HAVE_PARAFAFT.");
      return std::shared_ptr<FFTPlanInterface<double, NDim>>();
#endif
    }

  private:
#if defined(HAVE_MPI) && defined(HAVE_PARAFAFT)
    /** @brief Shared plan-creation path for any supported precision T. */
    template <typename T>
    std::shared_ptr<ParafaftPlanHolder<T, NDim>> make_plans(const MPICartesianGroup &group,
                                                            const FFTLayoutStruct<NDim> &layout)
    {
      auto globalSizes = layout.configurationSpace.getGlobalSizes();
      int globalShape[NDim];
      for (size_t i = 0; i < NDim; ++i)
        globalShape[i] = static_cast<int>(globalSizes[i]);

      auto parafaftObj =
          std::make_shared<parafaft::ParaFaFT_R2C<NDim, ParaFaFT_Backend<T>>>(globalShape, group.getBaseComm());

      return std::make_shared<ParafaftPlanHolder<T, NDim>>(group, parafaftObj);
    }
#endif
  };
} // namespace TempLat

#endif
