
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026
#include "TempLat/lattice/algebra/matrix3x3algebra/matrix3x3algebra.h"
#include "TempLat/lattice/algebra/complexalgebra/complexalgebra.h"

#include "TempLat/lattice/algebra/coordinates/wavenumber.h"
#include "TempLat/lattice/algebra/operators/operators.h"

#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/ndloop.h"


namespace TempLat
{

  struct AllMatrixMultiplyTester {
    static void Test(TDDAssertion &tdd);
  };


  void AllMatrixMultiplyTester::Test(TDDAssertion &tdd)
  {

    auto symTracelessStruct = ConstructSymTraceless(1., 2., 3., 4., 5., -5.);
    auto symStruct = ConstructSym(2., 4., 6., 8., 10., 12.);
    auto hermStruct = ConstructHerm(2., Complexify(3.,4.), Complexify(4.,5.), 8., Complexify(1.,2.), 12.);

    auto symsymtracelesstest = symStruct * symTracelessStruct;
    auto scalarsymtracelesstest = 4. * symTracelessStruct;
    auto hermsymtracelesstest = hermStruct * symTracelessStruct;

    /* Default is to fail: to remind yourself to implement something here. */
    tdd.verify(AlmostEqual(symsymtracelesstest.MatrixGet(1_c,1_c), 28));
    tdd.verify(AlmostEqual(symsymtracelesstest.MatrixGet(1_c,2_c), 50));
    tdd.verify(AlmostEqual(symsymtracelesstest.MatrixGet(1_c,3_c), -4));
    tdd.verify(AlmostEqual(symsymtracelesstest.MatrixGet(2_c,1_c), 50));
    tdd.verify(AlmostEqual(symsymtracelesstest.MatrixGet(2_c,2_c), 90));
    tdd.verify(AlmostEqual(symsymtracelesstest.MatrixGet(2_c,3_c), 2));
    tdd.verify(AlmostEqual(symsymtracelesstest.MatrixGet(3_c,1_c), 62));
    tdd.verify(AlmostEqual(symsymtracelesstest.MatrixGet(3_c,2_c), 112));
    tdd.verify(AlmostEqual(symsymtracelesstest.MatrixGet(3_c,3_c), 8));

    tdd.verify(AlmostEqual(scalarsymtracelesstest.SymTracelessGet(1_c,1_c), 4));
    tdd.verify(AlmostEqual(scalarsymtracelesstest.SymTracelessGet(1_c,2_c), 8));
    tdd.verify(AlmostEqual(scalarsymtracelesstest.SymTracelessGet(1_c,3_c), 12));
    tdd.verify(AlmostEqual(scalarsymtracelesstest.SymTracelessGet(2_c,1_c), 8));
    tdd.verify(AlmostEqual(scalarsymtracelesstest.SymTracelessGet(2_c,2_c), 16));
    tdd.verify(AlmostEqual(scalarsymtracelesstest.SymTracelessGet(2_c,3_c), 20));
    tdd.verify(AlmostEqual(scalarsymtracelesstest.SymTracelessGet(3_c,1_c), 12));
    tdd.verify(AlmostEqual(scalarsymtracelesstest.SymTracelessGet(3_c,2_c), 20));
    tdd.verify(AlmostEqual(scalarsymtracelesstest.SymTracelessGet(3_c,3_c), -20));

    tdd.verify(AlmostEqual(Real(hermsymtracelesstest.MatrixGet(1_c,1_c)), 20));
    tdd.verify(AlmostEqual(Real(hermsymtracelesstest.MatrixGet(1_c,2_c)), 36));
    tdd.verify(AlmostEqual(Real(hermsymtracelesstest.MatrixGet(1_c,3_c)), 1));
    tdd.verify(AlmostEqual(Real(hermsymtracelesstest.MatrixGet(2_c,1_c)), 22));
    tdd.verify(AlmostEqual(Real(hermsymtracelesstest.MatrixGet(2_c,2_c)), 43));
    tdd.verify(AlmostEqual(Real(hermsymtracelesstest.MatrixGet(2_c,3_c)), 44));
    tdd.verify(AlmostEqual(Real(hermsymtracelesstest.MatrixGet(3_c,1_c)), 42));
    tdd.verify(AlmostEqual(Real(hermsymtracelesstest.MatrixGet(3_c,2_c)), 72));
    tdd.verify(AlmostEqual(Real(hermsymtracelesstest.MatrixGet(3_c,3_c)), -43));

    tdd.verify(AlmostEqual(Imag(hermsymtracelesstest.MatrixGet(1_c,1_c)), 23));
    tdd.verify(AlmostEqual(Imag(hermsymtracelesstest.MatrixGet(1_c,2_c)), 41));
    tdd.verify(AlmostEqual(Imag(hermsymtracelesstest.MatrixGet(1_c,3_c)), -5));
    tdd.verify(AlmostEqual(Imag(hermsymtracelesstest.MatrixGet(2_c,1_c)), 2));
    tdd.verify(AlmostEqual(Imag(hermsymtracelesstest.MatrixGet(2_c,2_c)), 2));
    tdd.verify(AlmostEqual(Imag(hermsymtracelesstest.MatrixGet(2_c,3_c)), -22));
    tdd.verify(AlmostEqual(Imag(hermsymtracelesstest.MatrixGet(3_c,1_c)), -9));
    tdd.verify(AlmostEqual(Imag(hermsymtracelesstest.MatrixGet(3_c,2_c)), -18));
    tdd.verify(AlmostEqual(Imag(hermsymtracelesstest.MatrixGet(3_c,3_c)), -25));

    auto trace1 = multiplyTrace(symsymtracelesstest,symsymtracelesstest);
    tdd.verify(AlmostEqual(trace1.Get(), 13900));

    auto trace2 = multiplyTrace(hermsymtracelesstest,hermsymtracelesstest);
    tdd.verify(AlmostEqual(Real(trace2.Get()), 9898));
    tdd.verify(AlmostEqual(Imag(trace2.Get()), 0));

    auto trace3 = multiplyTrace(hermStruct,symTracelessStruct);
    tdd.verify(AlmostEqual(Real(trace3.Get()), 20));
    tdd.verify(AlmostEqual(Imag(trace3.Get()), 0));

    auto trace4 = multiplyTrace(symStruct,symTracelessStruct);
    tdd.verify(AlmostEqual(trace4.Get(), 126));

    {

      // Test multiplication between a symmetric matrix (equivalent to the real projector) and a symTraceless field in Fourier space
      constexpr size_t NDim = 3;
      using T = double;
      ptrdiff_t nGrid = 16, nGhost = 2;
      auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);
      toolBox->setVerbose();

      WaveNumber<NDim> k(toolBox);

      SymTracelessField<T, NDim> fs("fs", toolBox);

      Field<T, NDim> k1("k1", toolBox);
      Field<T, NDim> k2("k2", toolBox);
      Field<T, NDim> k3("k3", toolBox);

      k1.inFourierSpace() = k(1_c);
      k2.inFourierSpace() = k(2_c);
      k3.inFourierSpace() = k(3_c);

      fs.inFourierSpace()(1_c,1_c) = 2. / 3. * k1.inFourierSpace() - 1. / 3. * k2.inFourierSpace() - 1. /3. * k3.inFourierSpace();
      fs.inFourierSpace()(1_c,2_c) = k2.inFourierSpace();
      fs.inFourierSpace()(1_c,3_c) = k3.inFourierSpace();
      fs.inFourierSpace()(2_c,2_c) = - 1. / 3. * k1.inFourierSpace() + 2. / 3. * k2.inFourierSpace() - 1. /3. * k3.inFourierSpace();
      fs.inFourierSpace()(2_c,3_c) = k3.inFourierSpace();

      auto symmat = ConstructSym(k1.inFourierSpace(), k1.inFourierSpace(), k1.inFourierSpace(), k2.inFourierSpace(), k2.inFourierSpace(), k3.inFourierSpace());
      // auto hermmat = ConstructHerm(k1, Complexify(0., k1), Complexify(0., k1), k2, Complexify(0., k2), k3);

      auto symfs = symmat*fs.inFourierSpace();
      // auto hermfs = hermmat*fs;

      Field<T, NDim> s11("s11", toolBox);
      Field<T, NDim> s12("s12", toolBox);
      Field<T, NDim> s13("s13", toolBox);
      Field<T, NDim> s21("s21", toolBox);
      Field<T, NDim> s22("s22", toolBox);
      Field<T, NDim> s23("s23", toolBox);
      Field<T, NDim> s31("s31", toolBox);
      Field<T, NDim> s32("s32", toolBox);
      Field<T, NDim> s33("s33", toolBox);

      s11.inFourierSpace() = symfs(1_c, 1_c);
      s12.inFourierSpace() = symfs(1_c, 2_c);
      s13.inFourierSpace() = symfs(1_c, 3_c);
      s21.inFourierSpace() = symfs(2_c, 1_c);
      s22.inFourierSpace() = symfs(2_c, 2_c);
      s23.inFourierSpace() = symfs(2_c, 3_c);
      s31.inFourierSpace() = symfs(3_c, 1_c);
      s32.inFourierSpace() = symfs(3_c, 2_c);
      s33.inFourierSpace() = symfs(3_c, 3_c);

      {
        auto view1 = k1.inFourierSpace().getLocalNDHostView();
        auto view2 = k2.inFourierSpace().getLocalNDHostView();
        auto view3 = k3.inFourierSpace().getLocalNDHostView();
        auto view11 = s11.inFourierSpace().getLocalNDHostView();
        auto view12 = s12.inFourierSpace().getLocalNDHostView();
        auto view13 = s13.inFourierSpace().getLocalNDHostView();
        auto view21 = s21.inFourierSpace().getLocalNDHostView();
        auto view22 = s22.inFourierSpace().getLocalNDHostView();
        auto view23 = s23.inFourierSpace().getLocalNDHostView();
        auto view31 = s31.inFourierSpace().getLocalNDHostView();
        auto view32 = s32.inFourierSpace().getLocalNDHostView();
        auto view33 = s33.inFourierSpace().getLocalNDHostView();

        bool all_correct = true;
        NDLoop<NDim>(view11, [&](const auto... idx) {
          bool this_correct = true;
          this_correct &= AlmostEqual(view11(idx...), 2. / 3. * view1(idx...) * (view1(idx...) + view2(idx...) +  view3(idx...)) );
          this_correct &= AlmostEqual(view12(idx...), view1(idx...) * (-1. / 3. * view1(idx...) + 5. / 3. * view2(idx...) +  2. / 3. * view3(idx...)) );
          if (!this_correct) {
            std::cout << "Matrix operation test failed at index (";
            ((std::cout << idx << ", "), ...);
            std::cout << ") got (" << view11(idx...) << ", "
            << view12(idx...) << ", "
            << 2. / 3. * view1(idx...) * (view1(idx...) + view2(idx...) +  view3(idx...)) << ", "
            << view1(idx...) * (-1. / 3. * view1(idx...) + 5. / 3. * view2(idx...) +  2. / 3. * view3(idx...)) << ")\n";
          }
          all_correct = all_correct && this_correct;
        });
        tdd.verify(all_correct);
      }


        // Test multiplication between an hermitian matrix (equivalent to the complex projector) and a symTraceless field in Fourier space

        k1.inFourierSpace() = k(1_c);
        k2.inFourierSpace() = k(2_c);
        k3.inFourierSpace() = k(3_c);

        fs.inFourierSpace()(1_c,1_c) = 2. / 3. * k1.inFourierSpace() - 1. / 3. * k2.inFourierSpace() - 1. /3. * k3.inFourierSpace();
        fs.inFourierSpace()(1_c,2_c) = k2.inFourierSpace();
        fs.inFourierSpace()(1_c,3_c) = k3.inFourierSpace();
        fs.inFourierSpace()(2_c,2_c) = - 1. / 3. * k1.inFourierSpace() + 2. / 3. * k2.inFourierSpace() - 1. /3. * k3.inFourierSpace();
        fs.inFourierSpace()(2_c,3_c) = k3.inFourierSpace();

        auto hermmat = ConstructHerm(1. * k(1_c), Complexify(0., 1. * k(1_c)), Complexify(0., 1. * k(1_c)), 1. * k(2_c), Complexify(0., 1. * k(2_c)), 1. * k(3_c));

        auto hermfs = hermmat*fs.inFourierSpace();

        Field<T, NDim> h11r("h11r", toolBox);
        Field<T, NDim> h11i("h11i", toolBox);
        Field<T, NDim> h12r("h12r", toolBox);
        Field<T, NDim> h12i("h12i", toolBox);
        Field<T, NDim> h13r("h13r", toolBox);
        Field<T, NDim> h13i("h13i", toolBox);
        Field<T, NDim> h21r("h21r", toolBox);
        Field<T, NDim> h21i("h21i", toolBox);
        Field<T, NDim> h22r("h22r", toolBox);
        Field<T, NDim> h22i("h22i", toolBox);
        Field<T, NDim> h23r("h23r", toolBox);
        Field<T, NDim> h23i("h23i", toolBox);
        Field<T, NDim> h31r("h31r", toolBox);
        Field<T, NDim> h31i("h31i", toolBox);
        Field<T, NDim> h32r("h32r", toolBox);
        Field<T, NDim> h32i("h32i", toolBox);
        Field<T, NDim> h33r("h33r", toolBox);
        Field<T, NDim> h33i("h33i", toolBox);

        h11r.inFourierSpace() = Real(hermfs(1_c, 1_c));
        h11i.inFourierSpace() = Imag(hermfs(1_c, 1_c));
        h12r.inFourierSpace() = Real(hermfs(1_c, 2_c));
        h12i.inFourierSpace() = Imag(hermfs(1_c, 2_c));
        h13r.inFourierSpace() = Real(hermfs(1_c, 3_c));
        h13i.inFourierSpace() = Imag(hermfs(1_c, 3_c));
        h21r.inFourierSpace() = Real(hermfs(2_c, 1_c));
        h21i.inFourierSpace() = Imag(hermfs(2_c, 1_c));
        h22r.inFourierSpace() = Real(hermfs(2_c, 2_c));
        h22i.inFourierSpace() = Imag(hermfs(2_c, 2_c));
        h23r.inFourierSpace() = Real(hermfs(2_c, 3_c));
        h23i.inFourierSpace() = Imag(hermfs(2_c, 3_c));
        h31r.inFourierSpace() = Real(hermfs(3_c, 1_c));
        h31i.inFourierSpace() = Imag(hermfs(3_c, 1_c));
        h32r.inFourierSpace() = Real(hermfs(3_c, 2_c));
        h32i.inFourierSpace() = Imag(hermfs(3_c, 2_c));
        h33r.inFourierSpace() = Real(hermfs(3_c, 3_c));
        h33i.inFourierSpace() = Imag(hermfs(3_c, 3_c));

        {
          auto view1 = k1.inFourierSpace().getLocalNDHostView();
          auto view2 = k2.inFourierSpace().getLocalNDHostView();
          auto view3 = k3.inFourierSpace().getLocalNDHostView();
          auto view11r = h11r.inFourierSpace().getLocalNDHostView();
          auto view11i = h11i.inFourierSpace().getLocalNDHostView();
          auto view12r = h12r.inFourierSpace().getLocalNDHostView();
          auto view12i = h12i.inFourierSpace().getLocalNDHostView();
          auto view13r = h13r.inFourierSpace().getLocalNDHostView();
          auto view13i = h13i.inFourierSpace().getLocalNDHostView();
          auto view21r = h21r.inFourierSpace().getLocalNDHostView();
          auto view21i = h21i.inFourierSpace().getLocalNDHostView();
          auto view22r = h22r.inFourierSpace().getLocalNDHostView();
          auto view22i = h22i.inFourierSpace().getLocalNDHostView();
          auto view23r = h23r.inFourierSpace().getLocalNDHostView();
          auto view23i = h23i.inFourierSpace().getLocalNDHostView();
          auto view31r = h31r.inFourierSpace().getLocalNDHostView();
          auto view31i = h31i.inFourierSpace().getLocalNDHostView();
          auto view32r = h32r.inFourierSpace().getLocalNDHostView();
          auto view32i = h32i.inFourierSpace().getLocalNDHostView();
          auto view33r = h33r.inFourierSpace().getLocalNDHostView();
          auto view33i = h33i.inFourierSpace().getLocalNDHostView();

          bool all_correct = true;
          NDLoop<NDim>(view11r, [&](const auto... idx) {
            bool this_correct = true;
            this_correct &= AlmostEqual(view11r(idx...), view1(idx...) * (2. / 3. * view1(idx...) - 1. / 3. * view2(idx...) -  1. / 3. * view3(idx...)) );
            this_correct &= AlmostEqual(view11i(idx...), view1(idx...) * (view2(idx...) + view3(idx...)) );
            this_correct &= AlmostEqual(view21r(idx...), view2(idx...) * view2(idx...) );
            this_correct &= AlmostEqual(view21i(idx...), - view1(idx...) * (2. / 3. * view1(idx...) - 1. / 3. * view2(idx...) -  1. / 3. * view3(idx...)) + view2(idx...) * view3(idx...) );
            this_correct &= AlmostEqual(view31r(idx...), view3(idx...) * view3(idx...) );
            this_correct &= AlmostEqual(view31i(idx...), - view1(idx...) * (2. / 3. * view1(idx...) - 1. / 3. * view2(idx...) -  1. / 3. * view3(idx...)) - view2(idx...) * view2(idx...) );
            if (!this_correct) {
              std::cout << "Matrix operation test failed at index (";
              ((std::cout << idx << ", "), ...);
              std::cout << ") got (" << view11r(idx...) <<  "+ i * " << view11i(idx...) << ", "
              << view21r(idx...) <<  "+ i * " << view21i(idx...) << ", "
              << view31r(idx...) <<  "+ i * " << view31i(idx...) << ", "
              << view1(idx...) * (2. / 3. * view1(idx...) - 1. / 3. * view2(idx...) -  1. / 3. * view3(idx...)) <<  "+ i * " << view1(idx...) * (view2(idx...) + view3(idx...)) << ", "
              << view2(idx...) * view2(idx...) <<  "+ i * " << - view1(idx...) * (2. / 3. * view1(idx...) - 1. / 3. * view2(idx...) -  1. / 3. * view3(idx...)) + view2(idx...) * view3(idx...) <<  ", "
              << view3(idx...) * view3(idx...) <<  "+ i * " << - view1(idx...) * (2. / 3. * view1(idx...) - 1. / 3. * view2(idx...) -  1. / 3. * view3(idx...)) - view2(idx...) * view2(idx...) <<  ")\n";
            }
            all_correct = all_correct && this_correct;
          });
          tdd.verify(all_correct);
        }

    }

  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::AllMatrixMultiplyTester> test;
}
