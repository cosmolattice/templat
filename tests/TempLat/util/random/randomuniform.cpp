
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/util/random/randomuniform.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/almostequal.h"
#include <iomanip> // setprecision
#include "TempLat/parallel/device_iteration.h"

namespace TempLat
{
  struct RandomUniformTester {
    static void Test(TDDAssertion &tdd);
  };

  void RandomUniformTester::Test(TDDAssertion &tdd)
  {
    using T = double;
    {
      constexpr size_t N = 1000000;

      RandomUniform<T> prng("Hello CosmoLattice world!");

      double x = 0;

      // Let's draw two generations
      device::iteration::reduce<1>(
          "RandomUniformTester", {0}, {N},
          DEVICE_LAMBDA(device::IdxArray<1> i, double &sum) { sum += prng.get(i[0], i[0], 0); }, x);
      tdd.verify(AlmostEqual(x, 500167.407740739873));
      device::iteration::reduce<1>(
          "RandomUniformTester", {0}, {N},
          DEVICE_LAMBDA(device::IdxArray<1> i, double &sum) { sum += prng.get(i[0], i[0], 1); }, x);
      tdd.verify(AlmostEqual(x, 499917.024724008515));

      // And let's do that again:
      device::iteration::reduce<1>(
          "RandomUniformTester", {0}, {N},
          DEVICE_LAMBDA(device::IdxArray<1> i, double &sum) { sum += prng.get(i[0], i[0], 0); }, x);
      tdd.verify(AlmostEqual(x, 500167.407740739873));
      device::iteration::reduce<1>(
          "RandomUniformTester", {0}, {N},
          DEVICE_LAMBDA(device::IdxArray<1> i, double &sum) { sum += prng.get(i[0], i[0], 1); }, x);
      tdd.verify(AlmostEqual(x, 499917.024724008515));

      // Just a brief check, that all generated numbers are different
      device::memory::NDView<double, 1> a("a", 10);
      device::iteration::foreach<1>(
          "RandomUniformTester", {0}, {10},
          DEVICE_LAMBDA(device::IdxArray<1> i) { a(i[0]) = prng.get(i[0], i[0], 0); });
      auto a_host = device::memory::createMirrorView(a);
      device::memory::copyDeviceToHost(a, a_host.data());
      for (size_t i = 0; i < 9; ++i)
        tdd.verify(!AlmostEqual(a_host(i), a_host(i + 1)));

      // If we use the same state, all values should be identical
      device::memory::NDView<double, 1> b("b", 2);
      device::iteration::foreach<1>(
          "RandomUniformTester", {0}, {2}, DEVICE_LAMBDA(device::IdxArray<1> i) { b(i[0]) = prng.get(0, 0, 0); });
      auto b_host = device::memory::createMirrorView(b);
      device::memory::copyDeviceToHost(b, b_host.data());
      tdd.verify(AlmostEqual(b_host(0), b_host(1)));
    }
    {
      // Test saveState/loadState round-trip
      RandomUniform<T> rng("serialization_test");
      std::string savedState = rng.saveState();

      // Generate 1000 values after saving state
      std::vector<double> seq1;
      for (int i = 0; i < 1000; ++i) {
        seq1.push_back(rng.get(i, 0, 0));
      }

      // Restore state and generate again
      rng.loadState(savedState);
      std::vector<double> seq2;
      for (int i = 0; i < 1000; ++i) {
        seq2.push_back(rng.get(i, 0, 0));
      }

      // Verify sequences are identical
      bool sequencesMatch = true;
      for (int i = 0; i < 1000; ++i) {
        if (seq1[i] != seq2[i]) {
          sequencesMatch = false;
          break;
        }
      }
      tdd.verify(sequencesMatch, "saveState/loadState round-trip produces identical sequence");
    }
    {
      // --- Regression test for G6: a seed containing whitespace must survive save/loadState. ---
      // loadState reads the seed with `iss >> *mStringSeed`, which stops at the first whitespace, while
      // saveState writes the full string. A multi-word seed therefore round-trips to a different generator.
      RandomUniform<T> rng("Hello world"); // note the space
      RandomUniform<T> rng2("placeholder seed");
      rng2.loadState(rng.saveState());
      tdd.verify(rng2.getSeedString() == std::string("Hello world"),
                 "whitespace-containing seed survives saveState/loadState round-trip");
      tdd.verify(rng == rng2, "restored generator matches the original for a whitespace-containing seed");
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::RandomUniformTester> test;
}
