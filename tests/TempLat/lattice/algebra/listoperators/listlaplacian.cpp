/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026
#include "TempLat/lattice/algebra/listoperators/listoperators.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/field/collections/fieldcollection.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/operators/add.h"
#include "TempLat/lattice/algebra/operators/multiply.h"
#include "TempLat/lattice/algebra/operators/subtract.h"
#include "TempLat/lattice/measuringtools/averager.h"

namespace TempLat
{
  /** @brief The Laplacian of a field collection acts on each component as the scalar Laplacian does. */
  struct ListLaplacianTester {
    static void Test(TDDAssertion &tdd);
  };

  inline void ListLaplacianTester::Test(TDDAssertion &tdd)
  {
    constexpr size_t NDim = 3;
    auto toolBox = MemoryToolBox<NDim>::makeShared(8, 2);
    SpatialCoordinate<NDim> x(toolBox);

    FieldCollection<Field<double, NDim>, 2> fc("fc", toolBox);
    FieldCollection<Field<double, NDim>, 2> gc("gc", toolBox);
    fc(0_c) = x(1_c) * x(1_c) + 10. * x(2_c);
    fc(1_c) = 2. * x(1_c) * x(2_c) * x(3_c);

    auto deviation = [](const auto &a, const auto &b) { return average((a - b) * (a - b)); };

    gc = LatLapl(fc);
    tdd.verify(deviation(gc(0_c), LatLapl(fc(0_c))) == 0);
    tdd.verify(deviation(gc(1_c), LatLapl(fc(1_c))) == 0);

    gc = LatLapl(LatLapl(fc));
    tdd.verify(deviation(gc(1_c), LatLapl(LatLapl(fc(1_c)))) == 0);

    gc = fc + 0.5 * LatLapl(fc);
    tdd.verify(deviation(gc(0_c), fc(0_c) + 0.5 * LatLapl(fc(0_c))) == 0);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::ListLaplacianTester> test;
}
