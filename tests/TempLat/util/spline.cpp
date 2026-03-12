/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2025

#include "TempLat/util/spline.h"
#include "TempLat/parallel/device_iteration.h"
#include "TempLat/util/tdd/tdd.h"

#include <cmath>
#include <vector>

namespace TempLat
{

  template <typename T> struct SplineTester {
    static void Test(TDDAssertion &tdd);
  };

  // Evaluate a SplineData on device at a given point, returning the result to host.
  template <typename T> T evalOnDevice(SplineData<T> sd, T val)
  {
    T result{};
    device::iteration::reduce<1>(
        "SplineEval", device::IdxArray<1>{0}, device::IdxArray<1>{1},
        DEVICE_LAMBDA(const device::IdxArray<1> &, T &update) { update = sd(val); }, result);
    device::iteration::fence();
    return result;
  }

  // Evaluate a SplineData derivative on device, returning the result to host.
  template <typename T> T derivOnDevice(SplineData<T> sd, int order, T val)
  {
    T result{};
    device::iteration::reduce<1>(
        "SplineDeriv", device::IdxArray<1>{0}, device::IdxArray<1>{1},
        DEVICE_LAMBDA(const device::IdxArray<1> &, T &update) { update = sd.deriv(order, val); }, result);
    device::iteration::fence();
    return result;
  }

  template <typename T> void SplineTester<T>::Test(TDDAssertion &tdd)
  {
    const T tol = std::is_same_v<T, float> ? T(1e-4) : T(1e-10);

    // --- CPU fitting and host evaluation ---

    // Create data: y = x^2 on [0, 4]
    std::vector<T> xs = {T(0), T(1), T(2), T(3), T(4)};
    std::vector<T> ys = {T(0), T(1), T(4), T(9), T(16)};

    Spline<T> s(xs, ys, Spline<T>::cspline);

    // eval_host should pass through data points
    tdd.verify(std::fabs(s.eval_host(T(0)) - T(0)) < tol, "Host: passes through (0,0)");
    tdd.verify(std::fabs(s.eval_host(T(1)) - T(1)) < tol, "Host: passes through (1,1)");
    tdd.verify(std::fabs(s.eval_host(T(2)) - T(4)) < tol, "Host: passes through (2,4)");
    tdd.verify(std::fabs(s.eval_host(T(3)) - T(9)) < tol, "Host: passes through (3,9)");
    tdd.verify(std::fabs(s.eval_host(T(4)) - T(16)) < tol, "Host: passes through (4,16)");

    // Host interpolation at midpoints should be close to x^2
    tdd.verify(std::fabs(s.eval_host(T(0.5)) - T(0.25)) < T(0.1), "Host: interpolation at 0.5");
    tdd.verify(std::fabs(s.eval_host(T(1.5)) - T(2.25)) < T(0.1), "Host: interpolation at 1.5");
    tdd.verify(std::fabs(s.eval_host(T(2.5)) - T(6.25)) < T(0.1), "Host: interpolation at 2.5");

    // First derivative of x^2 is 2x
    tdd.verify(std::fabs(s.deriv_host(1, T(2)) - T(4)) < T(0.5), "Host: derivative at x=2");

    // --- Test linear spline ---
    Spline<T> slin(xs, ys, Spline<T>::linear);
    tdd.verify(std::fabs(slin.eval_host(T(2)) - T(4)) < tol, "Linear host: passes through (2,4)");
    tdd.verify(std::fabs(slin.eval_host(T(1.5)) - T(2.5)) < tol, "Linear host: midpoint interpolation");

    // --- Test hermite spline ---
    Spline<T> sherm(xs, ys, Spline<T>::cspline_hermite);
    tdd.verify(std::fabs(sherm.eval_host(T(2)) - T(4)) < tol, "Hermite host: passes through (2,4)");

    // --- Test solve ---
    auto roots = s.solve(T(4));
    bool found_2 = false;
    for (auto r : roots) {
      if (std::fabs(r - T(2)) < T(0.01)) found_2 = true;
    }
    tdd.verify(found_2, "solve(4) finds x=2");

    // --- Test make_monotonic ---
    {
      std::vector<T> mx = {T(0), T(1), T(2), T(3), T(4)};
      std::vector<T> my = {T(0), T(1), T(2), T(3), T(4)};
      Spline<T> sm(mx, my, Spline<T>::cspline, true);
      tdd.verify(sm.eval_host(T(0.5)) < sm.eval_host(T(1.5)), "Monotonic host: f(0.5) < f(1.5)");
      tdd.verify(sm.eval_host(T(1.5)) < sm.eval_host(T(2.5)), "Monotonic host: f(1.5) < f(2.5)");
    }

    // --- Device evaluation via SplineData ---
    // SplineData is captured by value into GPU lambdas and evaluated on-device.
    SplineData<T> sd = s.data();

    // Device eval should match eval_host
    tdd.verify(std::fabs(evalOnDevice(sd, T(0)) - s.eval_host(T(0))) < tol, "Device matches host at 0");
    tdd.verify(std::fabs(evalOnDevice(sd, T(1)) - s.eval_host(T(1))) < tol, "Device matches host at 1");
    tdd.verify(std::fabs(evalOnDevice(sd, T(2)) - s.eval_host(T(2))) < tol, "Device matches host at 2");
    tdd.verify(std::fabs(evalOnDevice(sd, T(1.5)) - s.eval_host(T(1.5))) < tol, "Device matches host at 1.5");
    tdd.verify(std::fabs(evalOnDevice(sd, T(3.7)) - s.eval_host(T(3.7))) < tol, "Device matches host at 3.7");

    // Device extrapolation
    tdd.verify(std::fabs(evalOnDevice(sd, T(-0.5)) - s.eval_host(T(-0.5))) < tol, "Device extrapolation left");
    tdd.verify(std::fabs(evalOnDevice(sd, T(5)) - s.eval_host(T(5))) < tol, "Device extrapolation right");

    // Device derivative should match host derivative
    tdd.verify(std::fabs(derivOnDevice(sd, 1, T(2)) - s.deriv_host(1, T(2))) < tol,
               "Device deriv(1) matches host at 2");
    tdd.verify(std::fabs(derivOnDevice(sd, 2, T(1.5)) - s.deriv_host(2, T(1.5))) < tol,
               "Device deriv(2) matches host at 1.5");
    tdd.verify(std::fabs(derivOnDevice(sd, 3, T(1.5)) - s.deriv_host(3, T(1.5))) < tol,
               "Device deriv(3) matches host at 1.5");
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SplineTester<double>> testDouble;
  TempLat::TDDContainer<TempLat::SplineTester<float>> testFloat;
} // namespace
