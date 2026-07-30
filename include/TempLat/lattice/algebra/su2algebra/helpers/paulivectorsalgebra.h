#ifndef TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_HELPERS_PAULIVECTORSALGEBRA_H
#define TEMPLAT_LATTICE_ALGEBRA_SU2ALGEBRA_HELPERS_PAULIVECTORSALGEBRA_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2025

#include "TempLat/parallel/device.h"

namespace TempLat
{

  /** @brief A class which implements the SU2 algebra at the single element level.
   *
   *
   * Unit test: ctest -R test-paulivectorsalgebra
   **/
  class PauliVectorsAlgebra
  {
  public:
    /* Put public methods here. These should change very little over time. */

    /**
     * @brief Multiplies two SU(2) elements
     *
     * @param res The result of the multiplication, which is modified in place. Must be an array of size 4.
     * @param cL The left element of the multiplication. Must be an array of size 4.
     * @param cR The right element of the multiplication. Must be an array of size 4.
     */
    template <typename ResArray, typename ArrayL, typename ArrayR>
      requires requires(ResArray r, ArrayL a, ArrayR b) {
        r[0];
        r[1];
        r[2];
        r[3];
        a[0];
        a[1];
        a[2];
        a[3];
        b[0];
        b[1];
        b[2];
        b[3];
      }
    DEVICE_INLINE_FUNCTION static void multiply_inplace(ResArray &res, const ArrayL &cL, const ArrayR &cR)
    {
      res[0] = cL[0] * cR[0] - cL[1] * cR[1] - cL[2] * cR[2] - cL[3] * cR[3];
      res[1] = cL[0] * cR[1] + cL[1] * cR[0] + cL[3] * cR[2] - cL[2] * cR[3];
      res[2] = cL[0] * cR[2] + cL[2] * cR[0] + cL[1] * cR[3] - cL[3] * cR[1];
      res[3] = cL[0] * cR[3] + cL[3] * cR[0] + cL[2] * cR[1] - cL[1] * cR[2];
    }

    /**
     * @brief Computes the exponential map for a given SU(2) algebra element.
     *
     * @param res The result array, which is modified in place. Must be an array of size 4.
     * @param alg The algebra array. Must be an array of size 4.
     */
    template <typename ResArray, typename AlgArray>
      requires requires(ResArray r, AlgArray a) {
        r[0];
        r[1];
        r[2];
        r[3];
        a[0];
        a[1];
        a[2];
      }
    DEVICE_INLINE_FUNCTION static void expmap_inplace(ResArray &res, const AlgArray &alg)
    {
      const auto a = device::sqrt(alg[0] * alg[0] + alg[1] * alg[1] + alg[2] * alg[2]);
      res[0] = device::cos(a);
      const auto sina = device::sin(a);
      // Guard: sin(a)/a → 1 as a → 0. Use direct comparison (device-compatible).
      if (a > decltype(a)(1e-15)) {
        const auto ratio = sina / a;
        res[1] = alg[0] * ratio;
        res[2] = alg[1] * ratio;
        res[3] = alg[2] * ratio;
      } else {
        res[1] = alg[0];
        res[2] = alg[1];
        res[3] = alg[2];
      }
    }

    /**
     * @brief Computes the exponential map for a given SU(2) algebra element, working in-place on the result array.
     *
     * @param res The result array, which is modified in place. Must be an array of size 4.
     * @param alg The algebra array. Must be an array of size 4.
     */
    template <typename ResArray>
      requires requires(ResArray r) {
        r[0];
        r[1];
        r[2];
        r[3];
      }
    DEVICE_INLINE_FUNCTION static void expmap_inplace(ResArray &res)
    {
      const auto a = device::sqrt(res[1] * res[1] + res[2] * res[2] + res[3] * res[3]);
      res[0] = device::cos(a);
      const auto sina = device::sin(a);
      // Guard: sin(a)/a → 1 as a → 0. Use direct comparison (device-compatible).
      if (a > decltype(a)(1e-15)) {
        const auto ratio = sina / a;
        res[1] = res[1] * ratio;
        res[2] = res[2] * ratio;
        res[3] = res[3] * ratio;
      } else {
        res[1] = res[1];
        res[2] = res[2];
        res[3] = res[3];
      }
    }
  };
} // namespace TempLat

#endif
