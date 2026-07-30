
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/lattice/algebra/helpers/getstring.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct GetStringTester {
    static void Test(TDDAssertion &tdd);
  };

  void GetStringTester::Test(TDDAssertion &tdd)
  {

    tdd.verify(std::is_same<std::string, decltype(GetString::get(complex<double>(1, 1)))>::value);

    tdd.verify(std::is_same<std::string, decltype(GetString::get((double)1))>::value);

    // tdd.verify( std::is_same<std::string, decltype(GetString::get(complex<double>(1, 1),2))>::value );
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::GetStringTester> test;
}
