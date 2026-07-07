/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2019
#include "TempLat/lattice/algebra/complexalgebra/arg.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfield.h"
#include "TempLat/util/ndloop.h"

namespace TempLat
{

  struct ImagTester {
    static void Test(TDDAssertion &tdd);
  };

  void ImagTester::Test(TDDAssertion &tdd)
  {
    struct MyStruct {
      double ComplexFieldGet(Tag<0> t) { return 3; };
      double ComplexFieldGet(Tag<1> t) { return 4; };
    };

    struct MyStruct {
      double ComplexFieldGet(Tag<0> t) { return 3; };
      double ComplexFieldGet(Tag<1> t) { return -4; };
    };

    MyStruct ms;
    MyStruct2 ms2;
    tdd.verify(AlmostEqual(arg(ms), 0.9272952180016122));
    tdd.verify(AlmostEqual(arg2(ms), 0.9272952180016122));

    tdd.verify(AlmostEqual(arg(ms), 5.355890089177974));
    tdd.verify(AlmostEqual(arg2(ms), -0.9272952180016122));

  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ImagTester> test;
}
