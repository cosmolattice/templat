#ifndef TEMPLAT_LATTICE_ALGEBRA_GAUGEALGEBRA_CENTEREDCOVARIANTDERIVATIVEO4_H
#define TEMPLAT_LATTICE_ALGEBRA_GAUGEALGEBRA_CENTEREDCOVARIANTDERIVATIVEO4_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020

// MATHEMATICS OF  COMPUTATION VOLUME 51,  NUMBER 184 OCTOBER 1988, PAGES 699-706 Generation of Finite Difference
// Formulas on Arbitrarily Spaced Grids By  Bengt Fornberg

#include "TempLat/util/assignabletuple.h"
#include "TempLat/lattice/algebra/listoperators/foldmultiply.h"
#include "TempLat/util/rangeiteration/make_list_tag.h"
#include "TempLat/lattice/algebra/helpers/getdx.h"

namespace TempLat
{
  /** @brief A class which computes a O(dx^4) discrete covariant derivative.
   *
   *
   * Unit test: ctest -R test-centeredcovariantderivativeo4
   *
   * @vocab-summary Centred gauge-covariant derivative on a five-point stencil, accurate to $O(dx^4)$.
   * @vocab-signature CenteredCovariantDerivativeO4<dim>(Us..., scalar)
   **/

  template <size_t dim, class... Args> auto CenteredCovariantDerivativeO4(Args... args)
  {
    constexpr size_t size = sizeof...(args);
    auto list = make_list(args...);

    auto scalar = list.getComp(Tag<size - 1>());

    auto dx = GetDx::getDx(scalar);

    auto UPlus = MakeVector(i, 1, dim, fold_multiply(MakeArray(j, 0, size - 2, list.getComp(j)(i))));
    auto UPlusPlus = MakeVector(i, 1, dim, fold_multiply(MakeArray(j, 0, size - 2, shift(list.getComp(j)(i), i))));
    auto UMinus =
        MakeVector(i, 1, dim, fold_multiply(MakeArray(j, 0, size - 2, dagger(shift(list.getComp(j)(i), -i)))));
    auto UMinusMinus = MakeVector(
        i, 1, dim, fold_multiply(MakeArray(j, 0, size - 2, dagger(shift(shift(list.getComp(j)(i), -i), -i)))));

    // The stencil weights are not exactly representable, so they are formed in the lattice's own
    // float type: as double literals they would promote the whole expression to fp64.
    using FT = std::decay_t<decltype(dx)>;
    auto Cov = MakeVector(i, 1, dim,
                          (1 / dx) * (-(FT(1) / FT(12)) * UPlusPlus(i) * shift(shift(scalar, i), i) +
                                      +(FT(2) / FT(3)) * UPlus(i) * shift(scalar, i) -
                                      (FT(2) / FT(3)) * UMinus(i) * shift(scalar, -i)) +
                              (FT(1) / FT(12)) * UMinusMinus(i) * shift(shift(scalar, -i), -i));
    return Cov;
  }
  /**
   * @vocab-summary The ungauged $O(dx^4)$ centred derivative.
   **/
  template <size_t dim, class T> auto CenteredDerivativeO4(T t) { return CenteredCovariantDerivativeO4<dim>(t); }
} // namespace TempLat

#endif
