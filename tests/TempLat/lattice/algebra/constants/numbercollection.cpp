
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2026

#include "TempLat/lattice/algebra/algebra.h"
#include "TempLat/lattice/algebra/constants/numbercollection.h"
#include "TempLat/util/almostequal.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/templatarray.h"

namespace TempLat
{

  struct NumberCollectionTester {
    static void Test(TDDAssertion &tdd);
  };

  void NumberCollectionTester::Test(TDDAssertion &tdd)
  {
    // --- Basic construction and Tag<I> access ---
    NumberCollection<double, 2> c;
    c(0_c).value = 1.0;
    c(1_c).value = 2.0;
    tdd.verify(c(0_c).value == 1.0);
    tdd.verify(c(1_c).value == 2.0);

    // --- getComp works (IsTempLatGettable) ---
    static_assert(IsTempLatGettable<0, NumberCollection<double, 2>>,
                  "NumberCollection must be IsTempLatGettable");
    tdd.verify(c.getComp(0_c).value == 1.0);
    tdd.verify(c.getComp(1_c).value == 2.0);

    // --- List-algebra: scalar * NumberCollection ---
    NumberCollection<double, 2> p;
    p(0_c).value = 3.0;
    p(1_c).value = 4.0;
    auto scaled = 2.0 * p;
    // scaled is a ListMultiplication<double, NumberCollection<double,2>>
    // getComp(Tag<0>) should give 2.0 * Number<double>{3.0} = Multiplication expr
    auto comp0 = GetComponent::get(scaled, 0_c);
    auto comp1 = GetComponent::get(scaled, 1_c);
    tdd.verify(AlmostEqual(DoEval::eval(comp0, size_t{0}), 6.0));
    tdd.verify(AlmostEqual(DoEval::eval(comp1, size_t{0}), 8.0));

    // --- Element-wise operator+= with another NumberCollection ---
    NumberCollection<double, 2> a;
    a(0_c).value = 1.0;
    a(1_c).value = 2.0;
    NumberCollection<double, 2> b;
    b(0_c).value = 10.0;
    b(1_c).value = 20.0;
    a += b;
    tdd.verify(AlmostEqual(a(0_c).value, 11.0));
    tdd.verify(AlmostEqual(a(1_c).value, 22.0));

    // --- Element-wise operator+= with scaled expression (evolver pattern) ---
    // fldS += coeff * piS
    NumberCollection<double, 2> fldS;
    fldS(0_c).value = 1.0;
    fldS(1_c).value = 2.0;
    NumberCollection<double, 2> piS;
    piS(0_c).value = 10.0;
    piS(1_c).value = 20.0;
    double coeff = 0.5;
    fldS += coeff * piS;
    tdd.verify(AlmostEqual(fldS(0_c).value, 6.0));
    tdd.verify(AlmostEqual(fldS(1_c).value, 12.0));

    // --- operator+=(ZeroType) is no-op ---
    NumberCollection<double, 2> z;
    z(0_c).value = 5.0;
    z(1_c).value = 7.0;
    z += ZeroType();
    tdd.verify(z(0_c).value == 5.0);
    tdd.verify(z(1_c).value == 7.0);

    // --- Assignment from TempLatArray / scalar (initialization pattern) ---
    // fldS0 / fStar where fldS0 is TempLatArray<T, NS>
    TempLatArray<double, 2> arr;
    arr[0] = 10.0;
    arr[1] = 20.0;
    double fStar = 5.0;
    NumberCollection<double, 2> nc;
    nc = arr / fStar;
    tdd.verify(AlmostEqual(nc(0_c).value, 2.0));
    tdd.verify(AlmostEqual(nc(1_c).value, 4.0));

    // --- Size 1 collection ---
    NumberCollection<double, 1> single;
    single(0_c).value = 42.0;
    tdd.verify(single(0_c).value == 42.0);
    single += 2.0 * single;
    tdd.verify(AlmostEqual(single(0_c).value, 126.0));

    // --- Copy assignment ---
    NumberCollection<double, 2> src;
    src(0_c).value = 100.0;
    src(1_c).value = 200.0;
    NumberCollection<double, 2> dst;
    dst = src;
    tdd.verify(dst(0_c).value == 100.0);
    tdd.verify(dst(1_c).value == 200.0);

    // --- float precision ---
    NumberCollection<float, 2> fc;
    fc(0_c).value = 1.5f;
    fc(1_c).value = 2.5f;
    fc += 0.5f * fc;
    tdd.verify(AlmostEqual(fc(0_c).value, 2.25f));
    tdd.verify(AlmostEqual(fc(1_c).value, 3.75f));
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::NumberCollectionTester> test;
}
