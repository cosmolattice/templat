
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/lattice/algebra/complexalgebra/complexfieldconjugate.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct ComplexFieldConjugateTester {
    static void Test(TDDAssertion &tdd);
  };

  void ComplexFieldConjugateTester::Test(TDDAssertion &tdd)
  {
    /* Default is to fail: to remind yourself to implement something here. */
    struct MyStruct {
      int ComplexFieldGet(Tag<0> t) const { return 1; }
      int ComplexFieldGet(Tag<1> t) const { return 2; }
    };
    tdd.verify(Real(conj(MyStruct())) == 1);
    tdd.verify(Imag(conj(MyStruct())) == -2);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ComplexFieldConjugateTester> test;
}
