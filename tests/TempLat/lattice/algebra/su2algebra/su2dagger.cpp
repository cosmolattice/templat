
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019
#include "TempLat/lattice/algebra/su2algebra/su2dagger.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/algebra/su2algebra/su2field.h"
#include "TempLat/lattice/algebra/su2algebra/su2multiply.h"

#include <cmath>
#include <algorithm>

namespace TempLat
{

  struct SU2DaggerTester {
    static void Test(TDDAssertion &tdd);
  };

  void SU2DaggerTester::Test(TDDAssertion &tdd)
  {
    auto toolBox = MemoryToolBox<3>::makeShared(16, 1);
    SU2Field<float, 3> su2_1("testSU2_1", toolBox);
    SU2Field<float, 3> result("resultSU2", toolBox);

    su2_1(1_c) = 1.0f;
    su2_1(2_c) = -1.5f;
    su2_1(3_c) = 2.0f;

    result = dagger(su2_1);

    auto f1view = result(1_c).getLocalNDHostView();
    auto f2view = result(2_c).getLocalNDHostView();
    auto f3view = result(3_c).getLocalNDHostView();

    bool all_true = true;
    for (size_t i = 0; i < f1view.extent(0); ++i)
      for (size_t j = 0; j < f1view.extent(1); ++j)
        for (size_t k = 0; k < f1view.extent(2); ++k) {
          all_true = all_true && (f1view(i, j, k) == -1.f);
          all_true = all_true && (f2view(i, j, k) == 1.5f);
          all_true = all_true && (f3view(i, j, k) == -2.f);

          if (!all_true) {
            std::cout << "Failed at index " << i << " " << j << " " << k << "\n";
            std::cout << "Values are: " << f1view(i, j, k) << " " << f2view(i, j, k) << " " << f3view(i, j, k) << "\n";
          }
        }
    // ---------------------------------------------------------------------------------------------
    // Equivalence of the dagger simplification rules (involution + cancelling anti-homomorphism) with
    // the definitional dagger. The reference materialises the product into a field and daggers that
    // leaf (which the rewrite rules do NOT touch), so any disagreement is a rule bug. The two sides do
    // a different multiply/add order, so compare within a float tolerance.
    // ---------------------------------------------------------------------------------------------
    SU2Field<float, 3> A("A", toolBox), B("B", toolBox);
    SU2Field<float, 3> prod("prod", toolBox), ref("ref", toolBox), tst("tst", toolBox);

    A(0_c) = 0.3f;  A(1_c) = 1.0f;  A(2_c) = -0.7f; A(3_c) = 0.5f;
    B(0_c) = -0.2f; B(1_c) = 0.4f;  B(2_c) = 1.1f;  B(3_c) = -0.9f;

    auto maxDiff = [](auto &x, auto &y) {
      float m = 0.f;
      auto cmp = [&](auto cx, auto cy) {
        auto vx = cx.getLocalNDHostView();
        auto vy = cy.getLocalNDHostView();
        for (size_t i = 0; i < vx.extent(0); ++i)
          for (size_t j = 0; j < vx.extent(1); ++j)
            for (size_t k = 0; k < vx.extent(2); ++k)
              m = std::max(m, std::abs(vx(i, j, k) - vy(i, j, k)));
      };
      cmp(x(0_c), y(0_c)); cmp(x(1_c), y(1_c)); cmp(x(2_c), y(2_c)); cmp(x(3_c), y(3_c));
      return m;
    };

    const float tol = 1e-4f;
    bool rules_ok = true;

    // (U^dagger)^dagger == U  (involution)
    tst = dagger(dagger(A));
    rules_ok = rules_ok && (maxDiff(tst, A) < tol);

    // dag(A . dag(B)) == (A . dag(B))^dagger   (right operand is a dagger -> cancels)
    prod = A * dagger(B);
    ref = dagger(prod);
    tst = dagger(A * dagger(B));
    rules_ok = rules_ok && (maxDiff(tst, ref) < tol);

    // dag(dag(A) . B) == (dag(A) . B)^dagger   (left operand is a dagger -> cancels)
    prod = dagger(A) * B;
    ref = dagger(prod);
    tst = dagger(dagger(A) * B);
    rules_ok = rules_ok && (maxDiff(tst, ref) < tol);

    // dag(dag(A) . dag(B)) == (dag(A) . dag(B))^dagger   (both operands are daggers)
    prod = dagger(A) * dagger(B);
    ref = dagger(prod);
    tst = dagger(dagger(A) * dagger(B));
    rules_ok = rules_ok && (maxDiff(tst, ref) < tol);

    // dag(A . B) == (A . B)^dagger   (no dagger operand -> rewrite must NOT fire, still correct)
    prod = A * B;
    ref = dagger(prod);
    tst = dagger(A * B);
    rules_ok = rules_ok && (maxDiff(tst, ref) < tol);

    if (!rules_ok)
      std::cout << "SU2 dagger simplification rules disagree with definitional dagger.\n";

    tdd.verify(all_true && rules_ok);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SU2DaggerTester> test;
}
