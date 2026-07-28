
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/lattice/algebra/conditional/conditionalunarygetter.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct ConditionalUnaryGetterTester {
    static void Test(TDDAssertion &tdd);
  };

  void ConditionalUnaryGetterTester::Test(TDDAssertion &tdd)
  {
    class MyClass
    {
    public:
      MyClass(int b) : a(b) {}

      DEVICE_INLINE_FUNCTION
      auto eval(device::Idx i) { return a; }

    private:
      int a;
    };

    tdd.verify(ConditionalUnaryGetter<MyClass> == true);
    tdd.verify(ConditionalUnaryGetter<double> == false);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ConditionalUnaryGetterTester> test;
}
