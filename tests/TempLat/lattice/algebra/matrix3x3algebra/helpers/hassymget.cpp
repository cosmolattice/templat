
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hassymget.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct HasSymGetTester {
    static void Test(TDDAssertion &tdd);
  };

  void HasSymGetTester::Test(TDDAssertion &tdd)
  {
    struct MyStruct {
      double SymGet(Tag<0> t) { return 87; };
    };
    struct MyStruct2 {
      double getComp(Tag<0> t) { return 87; };
    };
    struct MyStruct3 {
      double SymGet(Tag<1> t, Tag<1> t2) { return 64; };
    };

    tdd.verify(HasSymGet<MyStruct> == true);
    tdd.verify(HasSymGet<MyStruct2> == false);
    tdd.verify(HasSymGet<MyStruct3> == true);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::HasSymGetTester> test;
}
