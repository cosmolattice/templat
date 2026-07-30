
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026
#include "TempLat/lattice/algebra/matrix3x3algebra/matrix3x3algebra.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/almostequal.h"

namespace TempLat
{

  struct MatrixMatrixMultiplyTester {
    static void Test(TDDAssertion &tdd);
  };

  void MatrixMatrixMultiplyTester::Test(TDDAssertion &tdd)
  {
    // Two general (non-symmetric) 3x3 matrices, row-major.
    //   A = [[1,2,3],[4,5,6],[7,8,9]]   B = [[9,8,7],[6,5,4],[3,2,1]]
    auto A = ConstructMatrix3x3(1., 2., 3., 4., 5., 6., 7., 8., 9.);
    auto B = ConstructMatrix3x3(9., 8., 7., 6., 5., 4., 3., 2., 1.);

    auto C = A * B;

    // The symbolic MatrixGet(i,j) path is correct; use it as the reference for the flattened eval() path.
    // Hand-computed product (row-major, 0..8):
    const double expected[9] = {30, 24, 18, 84, 69, 54, 138, 114, 90};

    // --- Regression test for G1: the flattened device eval() path is wrong. ---
    // In MatrixMatrixMultiplication::eval(), result[8] is assigned three times while result[6] and result[7]
    // are never written, and the last column carries spurious minus signs. eval(0) therefore disagrees with
    // both the hand-computed product and the (correct) symbolic MatrixGet() values at indices 2, 5, 6, 7, 8.
    auto flat = C.eval(0);

    bool all_correct = true;
    for (int i = 0; i < 9; ++i) {
      if (!AlmostEqual(flat[i], expected[i])) {
        all_correct = false;
        say << "matA*matB eval() mismatch at flat index " << i << ": got " << flat[i] << ", expected "
            << expected[i] << "\n";
      }
    }
    tdd.verify(all_correct, "general 3x3 matA*matB flattened eval() matches the hand-computed product");
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::MatrixMatrixMultiplyTester> test;
}
