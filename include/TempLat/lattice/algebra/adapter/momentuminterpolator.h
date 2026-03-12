#ifndef TEMPLAT_LATTICE_ALGEBRA_ADAPTER_MOMENTUMINTERPOLATOR_H
#define TEMPLAT_LATTICE_ALGEBRA_ADAPTER_MOMENTUMINTERPOLATOR_H

#include <vector>
#include <algorithm>
#include <fstream>

#include "TempLat/lattice/algebra/coordinates/wavenumber.h"  // WaveNumber + FourierSite alias
#include "TempLat/lattice/memory/memorytoolbox.h"            // MemoryToolBox<NDim>
#include "TempLat/lattice/algebra/helpers/isvariadicindex.h" // IsVariadicNDIndex
#include "TempLat/util/spline.h"                              // spline

namespace TempLat
{
  template <typename T, size_t NDim> class MomentumInterpolator
  {
  public:
    using ToolBoxPtr = device::memory::host_ptr<MemoryToolBox<NDim>>;

    MomentumInterpolator(const std::vector<T> &kIn, const std::vector<T> &psIn, ToolBoxPtr toolBox, T kIRIn)
        : spline(kIn, psIn, Spline<T>::cspline), mToolBox(toolBox), ntilde(toolBox), ntilde_norm(ntilde.norm()),
          kIR(kIRIn), kMin(kIn.empty() ? T(0) : kIn.front()), kMax(kIn.empty() ? T(0) : kIn.back())
    {
    }

    // This is what TempLat operators (safeSqrt, +, *, etc.) want.
    template <typename... IDX>
      requires IsVariadicNDIndex<NDim, IDX...>
    DEVICE_FORCEINLINE_FUNCTION auto eval(const IDX &...idx) const
    {
      const auto kSite = ntilde_norm.eval(idx...) * kIR; // k = |ntilde| * kIR
      return spline(clampK(static_cast<T>(kSite)));
    }

    ToolBoxPtr getToolBox() const { return mToolBox; }

    std::string operatorString() const { return "Interpolator"; }

  private:
    Spline<T> spline;
    ToolBoxPtr mToolBox;

    FourierSite<NDim> ntilde;
    decltype(ntilde.norm()) ntilde_norm;
    T kIR, kMin, kMax;

    DEVICE_FORCEINLINE_FUNCTION
    T clampK(T k) const { return (k < kMin) ? kMin : (k > kMax) ? kMax : k; }
  };

} // namespace TempLat
#endif