/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026
#include "TempLat/lattice/measuringtools/accumulatortype.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  struct AccumulatorTypeTester {
    static void Test(TDDAssertion &tdd);
  };

  void AccumulatorTypeTester::Test(TDDAssertion &tdd)
  {
    /* Real and integral types accumulate in double, so that a lattice sum of N^3 single-precision
       terms does not lose its significance. */
    tdd.verify(std::is_same_v<AccumulatorType<float>, double>);
    tdd.verify(std::is_same_v<AccumulatorType<double>, double>);
    tdd.verify(std::is_same_v<AccumulatorType<int>, double>);

    /* Complex types promote component-wise: is_floating_point_v<complex<float>> is false, so a
       conditional on it alone would silently leave a Fourier-space sum in single precision. */
    tdd.verify(std::is_same_v<AccumulatorType<complex<float>>, complex<double>>);
    tdd.verify(std::is_same_v<AccumulatorType<complex<double>>, complex<double>>);

    /* bool counts as integral, so a count reduction accumulates in double too. */
    tdd.verify(std::is_same_v<AccumulatorType<bool>, double>);

    /* Anything non-arithmetic is left alone. */
    tdd.verify(std::is_same_v<AccumulatorType<std::string>, std::string>);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::AccumulatorTypeTester> test;
}
