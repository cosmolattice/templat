#ifndef TEMPLAT_UTIL_NUMERICALINTEGRATOR_H
#define TEMPLAT_UTIL_NUMERICALINTEGRATOR_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020

#include <vector>

namespace TempLat
{
  /** @brief A class which encapsulate some newto-cotes integration formula. Assume that the point are equispaced.
   *
   *
   * Unit test: ctest -R test-numericalintegrator
   **/
  class NumericalIntegrator
  {
  public:
    // Put public methods here. These should change very little over time.
    NumericalIntegrator() {}

    template <typename T> static T integrate(const std::vector<T> &vec, T dt)
    {
      const size_t n = vec.size();
      if (n < 2) return T(0); // no interval to integrate over

      const size_t m = n - 1;                                   // number of equispaced intervals
      const size_t simpsonIntervals = (m % 2 == 0) ? m : m - 1; // largest even interval count

      T res = T(0);
      // Composite Simpson over an even number of intervals (weights 1,4,2,4,...,4,1).
      if (simpsonIntervals >= 2) {
        T sum = vec[0] + vec[simpsonIntervals];
        for (size_t i = 1; i < simpsonIntervals; ++i)
          sum += (i % 2 == 1 ? T(4) : T(2)) * vec[i];
        res += dt / 3.0 * sum;
      }
      // If the interval count is odd, close the leftover interval with the trapezoidal rule.
      if (simpsonIntervals != m) res += dt / 2.0 * (vec[m - 1] + vec[m]);

      return res;
    }
  };

  template <typename T> T integrate(const std::vector<T> &vec, T dt) { return NumericalIntegrator::integrate(vec, dt); }

} // namespace TempLat

#endif
