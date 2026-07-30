
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/lattice/algebra/helpers/iscomplextype.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct IsComplexTypeTester {
    static void Test(TDDAssertion &tdd);
  };

  void IsComplexTypeTester::Test(TDDAssertion &tdd)
  {
    tdd.verify(IsComplexType<complex<float>>);
    tdd.verify(IsComplexType<complex<double>>);
    tdd.verify(IsComplexType<const complex<float>>);
    tdd.verify(IsComplexType<const complex<double>>);
    tdd.verify(IsComplexType<const complex<float> &>);
    tdd.verify(IsComplexType<const complex<double> &>);

    tdd.verify(!IsComplexType<int>);
    tdd.verify(!IsComplexType<std::array<int, 3>>);
    tdd.verify(!IsComplexType<float>);
    tdd.verify(!IsComplexType<double>);
    tdd.verify(!IsComplexType<float &>);
    tdd.verify(!IsComplexType<double &>);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::IsComplexTypeTester> test;
}
