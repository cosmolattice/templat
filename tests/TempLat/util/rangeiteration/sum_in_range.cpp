
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2019
#include "TempLat/util/rangeiteration/sum_in_range.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{
  struct this_should_compile {
    template <int n>
      requires(n > 5)
    static int func(Tag<n>)
    {
      return n;
    }
  };

  struct sum_in_range_Tester {
    static void Test(TDDAssertion &tdd);
  };

  void sum_in_range_Tester::Test(TDDAssertion &tdd)
  {
    /* Default is to fail: to remind yourself to implement something here. */
    tdd.verify(sum_in_range<1, 26>([](auto i) { return i; }) == (25 + 1) * 25 / 2);

    tdd.verify(sum_in_range<1, 1>([](auto i) { return this_should_compile::func(i); }).eval(0) == 0.);

    tdd.verify(Total(i, 1, 25, i) == (25 + 1) * 25 / 2);

    tdd.verify(Total(i, 1, 0, this_should_compile::func(i)).eval(0) == 0.);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::sum_in_range_Tester> test;
}
