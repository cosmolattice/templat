
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026
#include "TempLat/lattice/measuringtools/projectionhelpers/unbinnedradialprojectionresult.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/almostequal.h"

#include <tuple>

namespace TempLat
{

  struct UnbinnedRadialProjectionResultTester {
    static void Test(TDDAssertion &tdd);
  };

  void UnbinnedRadialProjectionResultTester::Test(TDDAssertion &tdd)
  {
    // Build a result with known bin positions and per-bin averages, bypassing finalize() (which needs MPI).
    UnbinnedRadialProjectionResult<double> r(4 /* nBins, unused by renormalizeBins */);

    auto datum = [](double avg) {
      RadialProjectionSingleDatum<double> d;
      d.average = avg;
      d.multiplicity = 1;
      return d;
    };

    // positions: 1, 4, 9, 16 ; averages: 10, 20, 30, 40
    r.push_back(std::make_tuple(1.0, datum(10.0)));
    r.push_back(std::make_tuple(4.0, datum(20.0)));
    r.push_back(std::make_tuple(9.0, datum(30.0)));
    r.push_back(std::make_tuple(16.0, datum(40.0)));

    r.renormalizeBins();

    // Expected per-bin widths (central differences; forward difference at the last bin):
    //   bin0: (p1 - p0)/2 - 0.5 = 1.0        -> 10 / 1 = 10
    //   bin1: (p2 - p0)/2       = 4.0        -> 20 / 4 = 5
    //   bin2: (p3 - p1)/2       = 6.0        -> 30 / 6 = 5
    //   bin3: (p3 - p2)         = 7.0        -> 40 / 7
    tdd.verify(AlmostEqual(std::get<1>(r[0]).average, 10.0), "first bin renormalized");
    tdd.verify(AlmostEqual(std::get<1>(r[1]).average, 5.0), "interior bin 1 renormalized");
    tdd.verify(AlmostEqual(std::get<1>(r[2]).average, 5.0), "interior bin 2 renormalized");

    // --- Regression test for A2: the last bin is left untouched (and index [size()] is read/written
    // out of bounds). The correct result divides the last bin by its width. ---
    tdd.verify(AlmostEqual(std::get<1>(r[3]).average, 40.0 / 7.0), "last bin renormalized (currently OOB / skipped)");

    // E3: toString() used to call getHeader()/toString() on the std::tuple element and did not compile if
    // instantiated. It should now render the stored (position, datum) tuples.
    tdd.verify(!r.toString().empty(), "toString renders the tuple entries");

    // Regression: integrate() reads binWidths, which only finalize() ever sized. On an instance whose
    // entries were filled in directly -- as above, to avoid needing a communicator -- that was an
    // out-of-bounds read. It should fall back to the unit widths finalize() would have installed.
    {
      UnbinnedRadialProjectionResult<double> u(2);
      u.push_back(std::make_tuple(1.0, datum(10.0)));
      u.push_back(std::make_tuple(2.0, datum(20.0)));

      // sum of average/position, each weighted by a unit bin width: 10/1 + 20/2.
      tdd.verify(AlmostEqual(u.integrate(false), 20.0), "integrate defaults to unit bin widths");
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::UnbinnedRadialProjectionResultTester> test;
}
