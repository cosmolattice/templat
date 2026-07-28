
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026
#include "TempLat/util/numericalintegrator.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/almostequal.h"

#include <vector>

namespace TempLat
{

  struct NumericalIntegratorTester {
    static void Test(TDDAssertion &tdd);
  };

  void NumericalIntegratorTester::Test(TDDAssertion &tdd)
  {
    // The exact 3-point Simpson rule works today: dt/3 * (1 + 4 + 1) = 2 for f == 1 over 2 intervals of dt = 1.
    std::vector<double> three{1.0, 1.0, 1.0};
    tdd.verify(AlmostEqual(NumericalIntegrator::integrate(three, 1.0), 2.0));

    // --- Regression test for G6: any sample count != 3 currently returns 0 silently. ---
    // The integral of the constant f == 1 over 5 equispaced samples (4 intervals of dt = 1) is 4, for any
    // correct quadrature rule. integrate() presently only handles size()==3 and returns 0 otherwise.
    std::vector<double> five{1.0, 1.0, 1.0, 1.0, 1.0};
    tdd.verify(AlmostEqual(NumericalIntegrator::integrate(five, 1.0), 4.0));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::NumericalIntegratorTester> test;
}
