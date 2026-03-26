
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026
#include "TempLat/lattice/algebra/matrix3x3algebra/matrix3x3algebra.h"
#include "TempLat/lattice/algebra/complexalgebra/complexalgebra.h"

#include "TempLat/lattice/algebra/coordinates/wavenumber.h"
#include "TempLat/lattice/algebra/operators/operators.h"

#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/ndloop.h"


namespace TempLat
{

  struct SymTracelessConjugateTester {
    static void Test(TDDAssertion &tdd);
  };


  void SymTracelessConjugateTester::Test(TDDAssertion &tdd)
  {

    using T = double;
    auto symTracelessStruct = ConstructSymTraceless(complex<T>(0., 1.), complex<T>(0., 2.), complex<T>(0., 3.), complex<T>(0., 4.), complex<T>(0., 5.), complex<T>(0., -5.));
    auto symTracelessStruct2 = ConstructSymTraceless(1., 2., 3., 4., 5., -5.);

    auto conjtest = conj(symTracelessStruct);
    auto conjtest2 = conj(symTracelessStruct2);

    /* Default is to fail: to remind yourself to implement something here. */
    tdd.verify(AlmostEqual(Real(conj(2)), 2));

    tdd.verify(AlmostEqual(Imag(conjtest(1_c,1_c)), -1));
    tdd.verify(AlmostEqual(Imag(conjtest(1_c,2_c)), -2));
    tdd.verify(AlmostEqual(Imag(conjtest(1_c,3_c)), -3));
    tdd.verify(AlmostEqual(Imag(conjtest(2_c,1_c)), -2));
    tdd.verify(AlmostEqual(Imag(conjtest(2_c,2_c)), -4));
    tdd.verify(AlmostEqual(Imag(conjtest(2_c,3_c)), -5));
    tdd.verify(AlmostEqual(Imag(conjtest(3_c,1_c)), -3));
    tdd.verify(AlmostEqual(Imag(conjtest(3_c,2_c)), -5));
    tdd.verify(AlmostEqual(Imag(conjtest(3_c,3_c)), 5));

    tdd.verify(AlmostEqual(Real(conjtest2(1_c,1_c)), 1.));
    tdd.verify(AlmostEqual(Real(conjtest2(1_c,2_c)), 2));
    tdd.verify(AlmostEqual(Real(conjtest2(1_c,3_c)), 3));
    tdd.verify(AlmostEqual(Real(conjtest2(2_c,1_c)), 2));
    tdd.verify(AlmostEqual(Real(conjtest2(2_c,2_c)), 4));
    tdd.verify(AlmostEqual(Real(conjtest2(2_c,3_c)), 5));
    tdd.verify(AlmostEqual(Real(conjtest2(3_c,1_c)), 3));
    tdd.verify(AlmostEqual(Real(conjtest2(3_c,2_c)), 5));
    tdd.verify(AlmostEqual(Real(conjtest2(3_c,3_c)), -5));

  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SymTracelessConjugateTester> test;
}
