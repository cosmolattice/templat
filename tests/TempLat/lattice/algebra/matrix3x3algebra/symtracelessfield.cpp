
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Jorge Baeza-Ballesteros,  Year: 2026

#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/matrix3x3algebra.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessfield.h"
#include "TempLat/util/ndloop.h"

namespace TempLat
{

  template <typename T, size_t NDim> struct SymTracelessTester {
    static void Test(TDDAssertion &tdd);
  };

  template <typename T, size_t NDim> inline void SymTracelessTester<T, NDim>::Test(TDDAssertion &tdd)
  {
    ptrdiff_t nGrid = 16, nGhost = 2;

    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);

    toolBox->setVerbose();

    // Test whether a transformation of the field forward and backward works.
    {
      SymTracelessField<T, NDim> original("original", toolBox);
      SpatialCoordinate x(toolBox);
      original.SymTracelessGet(Tag<0>()) = getVectorComponent(x, Tag<0>());
      original.SymTracelessGet(Tag<1>()) = getVectorComponent(x, Tag<1>());
      original.SymTracelessGet(Tag<2>()) = getVectorComponent(x, Tag<2>());
      original.SymTracelessGet(Tag<3>()) = getVectorComponent(x, Tag<0>()) + getVectorComponent(x, Tag<1>());
      original.SymTracelessGet(Tag<4>()) = getVectorComponent(x, Tag<0>()) + getVectorComponent(x, Tag<2>());
      original.updateGhosts();

      SymTracelessField<T, NDim> copy("copy", toolBox);
      copy = original;

      // force fourier transformation on copy
      copy.SymTracelessGet(Tag<0>()).getMemoryManager()->confirmFourierSpace();
      copy.SymTracelessGet(Tag<1>()).getMemoryManager()->confirmFourierSpace();
      copy.SymTracelessGet(Tag<2>()).getMemoryManager()->confirmFourierSpace();
      copy.SymTracelessGet(Tag<3>()).getMemoryManager()->confirmFourierSpace();
      copy.SymTracelessGet(Tag<4>()).getMemoryManager()->confirmFourierSpace();

      // force config transformation on copy
      copy.SymTracelessGet(Tag<0>()).getMemoryManager()->confirmConfigSpace();
      copy.SymTracelessGet(Tag<1>()).getMemoryManager()->confirmConfigSpace();
      copy.SymTracelessGet(Tag<2>()).getMemoryManager()->confirmConfigSpace();
      copy.SymTracelessGet(Tag<3>()).getMemoryManager()->confirmConfigSpace();
      copy.SymTracelessGet(Tag<4>()).getMemoryManager()->confirmConfigSpace();

      // update ghosts on copy
      copy.updateGhosts();

      auto original_host_0 = original.SymTracelessGet(Tag<0>()).getRawHostView();
      auto original_host_1 = original.SymTracelessGet(Tag<1>()).getRawHostView();
      auto original_host_2 = original.SymTracelessGet(Tag<2>()).getRawHostView();
      auto original_host_3 = original.SymTracelessGet(Tag<3>()).getRawHostView();
      auto original_host_4 = original.SymTracelessGet(Tag<4>()).getRawHostView();
      auto copy_host_0 = copy.SymTracelessGet(Tag<0>()).getRawHostView();
      auto copy_host_1 = copy.SymTracelessGet(Tag<1>()).getRawHostView();
      auto copy_host_2 = copy.SymTracelessGet(Tag<2>()).getRawHostView();
      auto copy_host_3 = copy.SymTracelessGet(Tag<3>()).getRawHostView();
      auto copy_host_4 = copy.SymTracelessGet(Tag<4>()).getRawHostView();

      bool backforthWorks = true;
      for (size_t i = 0; i < original_host_0.size(); ++i) {
        backforthWorks = backforthWorks && AlmostEqual(original_host_0[i], copy_host_0[i]);
        backforthWorks = backforthWorks && AlmostEqual(original_host_1[i], copy_host_1[i]);
        backforthWorks = backforthWorks && AlmostEqual(original_host_2[i], copy_host_2[i]);
        backforthWorks = backforthWorks && AlmostEqual(original_host_3[i], copy_host_3[i]);
        backforthWorks = backforthWorks && AlmostEqual(original_host_4[i], copy_host_4[i]);
      }
      tdd.verify(backforthWorks);
    }
    // ------------------------------------------------------------------------------------------
    {

      SymTracelessField<T, NDim> phi("phi", toolBox);
      SymTracelessField<T, NDim> chi("chi", toolBox);

      auto field_tester = [&](SymTracelessField<T, NDim> &f, auto op, std::array<double, 5> expected) {
        f = op;
        auto view0 = f.SymTracelessGet(0_c).getLocalNDHostView();
        auto view1 = f.SymTracelessGet(1_c).getLocalNDHostView();
        auto view2 = f.SymTracelessGet(2_c).getLocalNDHostView();
        auto view3 = f.SymTracelessGet(3_c).getLocalNDHostView();
        auto view4 = f.SymTracelessGet(4_c).getLocalNDHostView();

        bool all_correct = true;
        NDLoop<NDim>(view0, [&](const auto... idx) {
          bool this_correct = true;
          this_correct &= AlmostEqual(view0(idx...), expected[0]);
          this_correct &= AlmostEqual(view1(idx...), expected[1]);
          this_correct &= AlmostEqual(view2(idx...), expected[2]);
          this_correct &= AlmostEqual(view3(idx...), expected[3]);
          this_correct &= AlmostEqual(view4(idx...), expected[4]);
          if (!this_correct) {
            std::cout << "SymTracelessField operation test failed at index (";
            ((std::cout << idx << ", "), ...);
            std::cout << ") got (" << view0(idx...) << ", " << view1(idx...) << ", " << view2(idx...) << ", " << view3(idx...) << ", " << view4(idx...) << ") , expected (" << expected[0] << ", " << expected[1] << ", " << expected[2] << ", " << expected[3] << ", " << expected[4] << ")\n";
          }
          all_correct = all_correct && this_correct;
        });
        tdd.verify(all_correct);
      };

      getComponent(chi, 0_c) = 1;
      getComponent(chi, 2_c) = 2;
      getComponent(chi, 3_c) = 5;
      getComponent(phi, 1_c) = 8;
      getComponent(phi, 2_c) = 6;
      getComponent(phi, 4_c) = 3;

      field_tester(chi, chi, {1, 0, 2, 5, 0});
      field_tester(phi, chi + chi, {2, 0, 4, 10, 0});
      // field_tester(phi, chi * chi, {1, 0, 4, 25, 0});
      field_tester(phi, chi - chi, {0, 0, 0, 0, 0});

      /*
      field_tester(phi, cos(chi), cos(2));
      field_tester(phi, sin(chi), sin(2));
      field_tester(phi, tanh(chi), tanh(2));
      field_tester(phi, sqrt(chi), sqrt(2));
      field_tester(phi, log(chi), log(2));
      field_tester(phi, abs(chi), abs(2));
      field_tester(phi, asinh(chi), asinh(2));
    }

    // ------------------------------------------------------------------------------------------

    {
      Field<T, NDim> phi("phi", toolBox);
      Field<T, NDim> chi("chi", toolBox);
      Field<T, NDim> psi("psi", toolBox);

      std::cout << "Layout info: " << toolBox->mLayouts.getConfigSpaceLayout() << "\n";

      phi.inFourierSpace();
      tdd.verify(phi.mManager->isFourierSpace());

      WaveNumber k(toolBox);
      phi.inFourierSpace() = k.norm2() * RandomGaussianField<T, NDim>("Hoi", toolBox);

      // just manipulated phi(k), so it must still be in Fourier space, and ghosts are stale.
      tdd.verify(phi.mManager->isFourierSpace());
      tdd.verify(phi.mManager->areGhostsStale());

      chi = LatticeLaplacian<NDim, decltype(phi)>(phi);
      // just manipulated chi(x), so it must still be in configuration space, and ghosts are stale.
      tdd.verify(!chi.mManager->isFourierSpace());
      tdd.verify(chi.mManager->isConfigSpace());
      tdd.verify(chi.mManager->areGhostsStale());
      tdd.verify(!phi.mManager->isFourierSpace());
      tdd.verify(!phi.mManager->areGhostsStale());

      // neat consequence of the implementation: an expression actually evaluates to a specific type. Keeping that
      // instance, without passing it to an assignment operator, is simply the compiled expression. So we can do stuff
      // with it.
      auto potential = 0.5 * phi * phi + 42 * chi * chi * phi * phi - chi + (-chi);

      // Stuff we can do:
      say << "Potential2: " << potential.toString() << "\n";

      // and awesomer:
      // auto dVdPhi = potential.d(phi);

      // say << "dPotential/dphi: " << dVdPhi.toString() << "\n";
    */
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::SymTracelessTester<double, 3>> test;
}
