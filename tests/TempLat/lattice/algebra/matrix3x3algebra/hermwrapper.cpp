
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026
#include "TempLat/lattice/algebra/matrix3x3algebra/hermwrapper.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfieldconjugate.h"
#include "TempLat/lattice/algebra/complexalgebra/complexwrapper.h"
#include "TempLat/lattice/algebra/complexalgebra/real.h"
#include "TempLat/lattice/algebra/complexalgebra/imag.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct SymmetricTracelessWrapperTester {
    static void Test(TDDAssertion &tdd);
  };

  void SymmetricTracelessWrapperTester::Test(TDDAssertion &tdd)
  {

    auto test = ConstructHerm(1., Complexify(1, 2), Complexify(1, 3), 4., Complexify(1, 6), 6.);
    /* Default is to fail: to remind yourself to implement something here. */
    tdd.verify(AlmostEqual(test.HermGet(0_c), 1));
    tdd.verify(AlmostEqual(Real(test.HermGet(1_c)), 1));
    tdd.verify(AlmostEqual(Imag(test.HermGet(1_c)), 2));
    tdd.verify(AlmostEqual(Real(test.HermGet(2_c)), 1));
    tdd.verify(AlmostEqual(Imag(test.HermGet(2_c)), 3));
    tdd.verify(AlmostEqual(test.HermGet(3_c), 4));
    tdd.verify(AlmostEqual(Real(test.HermGet(4_c)), 1));
    tdd.verify(AlmostEqual(Imag(test.HermGet(4_c)), 6));
    tdd.verify(AlmostEqual(test.HermGet(5_c), 6));

    tdd.verify(AlmostEqual(test.HermGet(1_c,1_c), 1));
    tdd.verify(AlmostEqual(Real(test.HermGet(1_c,2_c)), 1));
    tdd.verify(AlmostEqual(Imag(test.HermGet(1_c,2_c)), 2));
    tdd.verify(AlmostEqual(Real(test.HermGet(1_c,3_c)), 1));
    tdd.verify(AlmostEqual(Imag(test.HermGet(1_c,3_c)), 3));
    tdd.verify(AlmostEqual(Real(test.HermGet(2_c,1_c)), 1));
    tdd.verify(AlmostEqual(Imag(test.HermGet(2_c,1_c)), -2));
    tdd.verify(AlmostEqual(test.HermGet(2_c,2_c), 4));
    tdd.verify(AlmostEqual(Real(test.HermGet(2_c,3_c)), 1));
    tdd.verify(AlmostEqual(Imag(test.HermGet(2_c,3_c)), 6));
    tdd.verify(AlmostEqual(Real(test.HermGet(3_c,1_c)), 1));
    tdd.verify(AlmostEqual(Imag(test.HermGet(3_c,1_c)), -3));
    tdd.verify(AlmostEqual(Real(test.HermGet(3_c,2_c)), 1));
    tdd.verify(AlmostEqual(Imag(test.HermGet(3_c,2_c)), -6));
    tdd.verify(AlmostEqual(test.HermGet(3_c,3_c), 6));

  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SymmetricTracelessWrapperTester> test;
}
