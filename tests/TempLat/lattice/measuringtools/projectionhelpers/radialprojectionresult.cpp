
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionresult.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/log/saycomplete.h"
#include "TempLat/util/powr.h"
#include "TempLat/parallel/device_iteration.h"

namespace TempLat
{

  template <typename T> struct RadialProjectionResultTester {
    static void Test(TDDAssertion &tdd);
  };

  template <typename T> inline void RadialProjectionResultTester<T>::Test(TDDAssertion &tdd)
  {
    /* (nBins, deltakBins) — the bin spacing only enters integrate(), which this test does not exercise. */
    RadialProjectionResult<T> one(10, 1), two(12, 1), three(10, 1);

    // tdd.verify(Throws<RadialProjectionResultSizeException>([&]() { one += two; }));

    /* dummy data */
    // for (device::Idx i = 0, iEnd = three.size(); i < iEnd; ++i) {
    //   three.add(i, 2 * i, 2 * i);
    // }

    int size = 1;
#ifdef HAVE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    device::iteration::foreach<1>(
        "RadialProjectionResultTest", {0}, {10}, DEVICE_LAMBDA(const device::IdxArray<1> i) {
          three.add_device(i[0], 2 * i[0], 2 * i[0]);
          three.add_device(i[0], 2 * i[0] + 0.5, 2 * i[0]);
          three.add_device(i[0], 2 * i[0] - 0.5, 2 * i[0]);
        });
    three.finalize(MPICommReference());
    tdd.verify(three.size() == 10);
    {
      bool allRight = true;
      for (device::Idx i = 0, iEnd = three.size(); i < iEnd; ++i) {
        const double expectedVariance =
            abs(powr<2>(2. * i) * 2. / 3. - powr<2>(2 * i + 0.5) / 3. - powr<2>(2 * i - 0.5) / 3.);
        allRight = allRight && AlmostEqual(three[i].getValue().average, 2 * i);
        allRight = allRight && AlmostEqual(three[i].getValue().sampleVariance, expectedVariance);
        allRight = allRight && three[i].getValue().multiplicity == 3 * size;
        if (!allRight) {
          say << "Bin " << i << ": "
              << "Average: " << three[i].getValue().average
              << ", Sample Variance: " << three[i].getValue().sampleVariance
              << ", Multiplicity: " << three[i].getValue().multiplicity << "\n"
              << "Expected: "
              << "Average: " << 2 * i << ", Sample Variance: " << expectedVariance << ", Multiplicity: " << 3 * size
              << "\n";

          break;
        }
      }
      tdd.verify(allRight);
    }
    /*
      one += three;
      one += three;

      bool allRight = true;
      tdd.verify(one.size() == three.size());

      for (device::Idx i = 0, iEnd = one.size(); i < iEnd; ++i) {
        allRight = allRight && one[i].getValue().average == 2 * three.mValues.mAverages[i];
        allRight = allRight && one[i].getValue().sampleVariance == 2 * three.mValues.mVariances[i];
        allRight = allRight && one[i].getValue().multiplicity == 2 * three.mMultiplicities[i];

        say << "Bin " << i << ": "
            << "Average: " << one[i].getValue().average << ", Sample Variance: " << one[i].getValue().sampleVariance
            << ", Multiplicity: " << one[i].getValue().multiplicity << "\n";
      }

      tdd.verify(allRight);
    */
    /* Regression: the converting constructor omitted mDeltakBins from its initializer list, leaving it
       indeterminate, while integrate() scales its total by it. Build the finalized shape by hand here --
       push_back into the parent vector and set the bin bounds through the public accessor -- so the check
       does not need a communicator. */
    {
      auto datum = [](double avg) {
        RadialProjectionSingleDatum<double> d;
        d.average = avg;
        d.multiplicity = 1;
        return d;
      };

      const double deltakBins = 2.5;
      RadialProjectionResult<double> src(2, deltakBins);
      src.getCentralBinBounds()[0] = 1.0;
      src.getCentralBinBounds()[1] = 2.0;
      src.push_back(RadialProjectionSingleBinAndValue<double>(datum(1.0), datum(10.0)));
      src.push_back(RadialProjectionSingleBinAndValue<double>(datum(2.0), datum(20.0)));

      /* integrate(true) sums average/centralBinBound over the bins and scales by deltakBins. */
      const double expected = (10.0 / 1.0 + 20.0 / 2.0) * deltakBins;
      tdd.verify(AlmostEqual(src.integrate(true), expected), "integrate scales by deltakBins");

      /* float != double, so this selects the converting constructor rather than the copy constructor. */
      RadialProjectionResult<float> converted(src);
      tdd.verify(AlmostEqual(converted.integrate(true), expected, 1e-5),
                 "converting constructor carries deltakBins across");
    }

    /* test that this compiles */
    for (auto &&it : one) {
      it.getValue().average *= 1;
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::RadialProjectionResultTester<double>> test;
}
