
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler,  Year: 2025
#include "TempLat/lattice/algebra/helpers/getndim.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/field/field.h"

namespace TempLat
{

  struct GetNDimTester {
    static void Test(TDDAssertion &tdd);
  };

  void GetNDimTester::Test(TDDAssertion &tdd)
  {
    tdd.verify(GetNDim::get<Field<double, 1>>() == 1);
    tdd.verify(GetNDim::get<Field<double, 2>>() == 2);
    tdd.verify(GetNDim::get<Field<double, 3>>() == 3);
    tdd.verify(GetNDim::get<Field<double, 4>>() == 4);
    tdd.verify(GetNDim::get<Field<double, 5>>() == 5);
    tdd.verify(GetNDim::get<Field<double, 6>>() == 6);

    tdd.verify(GetNDim::get<double>() == 0);
    tdd.verify(GetNDim::get<complex<double>>() == 0);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::GetNDimTester> test;
}
