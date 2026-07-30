
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020
#include "TempLat/lattice/algebra/operators/arg.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct ArgTester {
    static void Test(TDDAssertion &tdd);
  };

  void ArgTester::Test(TDDAssertion &tdd)
  {
    struct MyStructR {
      DEVICE_INLINE_FUNCTION
      double eval(int i) const { return (double)i / 4.0; }
    };

    struct MyStructC {
      MyStructR ComplexFieldGet(Tag<0> t) { return r; }
      MyStructR ComplexFieldGet(Tag<1> t) { return r; }

      MyStructR r;
    };

    auto t1 = MyStructC();
    auto argu = arg(t1.ComplexFieldGet(0_c), t1.ComplexFieldGet(1_c));

    tdd.verify(AlmostEqual(argu.eval(0), 0));
    tdd.verify(AlmostEqual(argu.eval(1), 0.7853981633974483));
    tdd.verify(AlmostEqual(argu.eval(-1), 3.926990816987241));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ArgTester> test;
}
