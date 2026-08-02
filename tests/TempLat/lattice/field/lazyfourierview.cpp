/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2026

/* Field holds its FourierView behind a lazily created host_ptr rather than by value, so that
 * a field which never enters Fourier space -- which is every field in the MC and cooling
 * paths -- pays neither its bytes nor its reference count on the copies that build every
 * expression tree.
 *
 * The existing Fourier tests (fieldviewfourier, field, randomgaussianfield, radialprojector)
 * all exercise Fourier-space *functionality*, and they passed before this change and after
 * it. None of them exercises the *laziness*: that the handle starts null, that it is
 * materialised exactly once, that copies share it or not depending on when they were taken,
 * and that destroying a null handle is safe. Those are the invariants that would break
 * silently, so they are what this file pins.
 */

#include "TempLat/lattice/field/field.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/util/ndloop.h"

#include <vector>

namespace TempLat
{

  template <typename T, size_t NDim> struct LazyFourierViewTester {
    static void Test(TDDAssertion &tdd);

    /** @brief Stands in for FieldCollection::operator(), which returns a Field BY VALUE
     *  (fieldcollection.h). Callers then write model.fldCS(i)(0_c).inFourierSpace() = ...,
     *  binding a reference into a temporary. */
    static Field<T, NDim> byValue(const Field<T, NDim> &f) { return f; }
  };

  template <typename T, size_t NDim> inline void LazyFourierViewTester<T, NDim>::Test(TDDAssertion &tdd)
  {
    /* The whole sharing argument rests on a FourierView being unable to be overwritten
     * through the reference inFourierSpace() hands out: the implicit copy assignment is
     * deleted by the const LatticeParameters member, and the operator=(R&&) template does
     * not win against it. If that ever stops holding, two Fields sharing a view stop being
     * equivalent to two Fields owning identical ones, and the change becomes observable. */
    static_assert(!std::is_copy_assignable_v<FourierView<T, NDim>>,
                  "FourierView must not be copy-assignable: Field copies may share one.");

    const device::Idx nGrid = 8, nGhost = 1;
    auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, nGhost);

    // ---- materialised exactly once, and stable across repeated calls -----------------
    {
      Field<T, NDim> a("lazy_once", toolBox);

      const FourierView<T, NDim> *first = &a.inFourierSpace();
      tdd.verify(first == &a.inFourierSpace());

      bool stable = true;
      for (int i = 0; i < 32; ++i)
        stable = stable && (&a.inFourierSpace() == first);
      tdd.verify(stable);
    }

    // ---- the handle really does start null -------------------------------------------
    /* Observed through its only visible consequence: a copy taken BEFORE the original is
     * ever used in Fourier space carries a null handle, so it materialises its own view;
     * a copy taken AFTER shares the original's. If the handle were eagerly created in the
     * constructor, both copies would share and the first check would fail. */
    {
      Field<T, NDim> a("lazy_null", toolBox);
      Field<T, NDim> early = a; // copied while both handles are null

      const FourierView<T, NDim> *va = &a.inFourierSpace();
      const FourierView<T, NDim> *ve = &early.inFourierSpace();
      tdd.verify(va != ve);

      Field<T, NDim> late = a; // copied after a's handle exists
      tdd.verify(&late.inFourierSpace() == va);
    }

    // ---- a field that never enters Fourier space survives being copied a lot ---------
    /* This is the hot path: expression-template nodes store operands by value, so a leaf is
     * copy-constructed once per level of tree depth and destroyed again. Every one of those
     * copies and destructions goes through the null branch of host_ptr, which returns early
     * instead of touching a reference count. A refcount bug there would corrupt or free the
     * field's memory, so assert the data is still intact afterwards. */
    {
      Field<T, NDim> a("lazy_neverfourier", toolBox);
      const T value = T(2.75);
      a = value;

      {
        std::vector<Field<T, NDim>> copies;
        for (int i = 0; i < 16; ++i)
          copies.push_back(a);
      } // all 16 destroyed here, every handle null

      tdd.verify(a.mManager->isConfigSpace());

      auto host = a.getLocalNDHostView();
      bool intact = true;
      NDLoop<NDim>(host, [&](const auto &...idx) { intact = intact && AlmostEqual(host(idx...), value); });
      tdd.verify(intact);
    }

    // ---- config -> Fourier -> config round trip --------------------------------------
    {
      Field<T, NDim> a("lazy_roundtrip", toolBox);
      const T value = T(3.25);
      a = value;

      a.inFourierSpace(); // forces R2C and materialises the view
      tdd.verify(a.mManager->isFourierSpace());

      a.mManager->confirmConfigSpace(); // C2R back
      tdd.verify(a.mManager->isConfigSpace());

      auto host = a.getLocalNDHostView();
      bool same = true;
      NDLoop<NDim>(host, [&](const auto &...idx) { same = same && AlmostEqual(host(idx...), value); });
      tdd.verify(same);
    }

    // ---- writing through a temporary Field lands in the original's memory ------------
    /* The FieldCollection::operator() pattern. Before this change the reference pointed at a
     * subobject of the temporary; now it points at a heap view the temporary solely owns.
     * Either way the write goes through the shared MemoryManager, so the result must be
     * indistinguishable from writing through a named field. */
    {
      Field<T, NDim> viaTemporary("lazy_tmp", toolBox);
      Field<T, NDim> viaNamed("lazy_named", toolBox);
      const T value = T(1.5);

      byValue(viaTemporary).inFourierSpace() = value; // reference into a temporary
      viaNamed.inFourierSpace() = value;

      auto hTmp = viaTemporary.inFourierSpace().getLocalNDHostView();
      auto hNamed = viaNamed.inFourierSpace().getLocalNDHostView();

      bool same = true;
      NDLoop<NDim>(hTmp, [&](const auto &...idx) { same = same && AlmostEqual(hTmp(idx...), hNamed(idx...)); });
      tdd.verify(same);
    }

    // ---- a shared view and an independently materialised one agree -------------------
    /* The sharing introduced by the handle is only safe if two views built from the same
     * Field are interchangeable. Write through one, read through the other. */
    {
      Field<T, NDim> a("lazy_shared", toolBox);
      Field<T, NDim> reference("lazy_shared_ref", toolBox);
      const T value = T(4.5);

      a.inFourierSpace() = value; // materialise, then share
      Field<T, NDim> sharer = a;
      reference.inFourierSpace() = value;

      auto hShared = sharer.inFourierSpace().getLocalNDHostView();
      auto hRef = reference.inFourierSpace().getLocalNDHostView();

      bool same = true;
      NDLoop<NDim>(hShared, [&](const auto &...idx) { same = same && AlmostEqual(hShared(idx...), hRef(idx...)); });
      tdd.verify(same);
    }
  }

} // namespace TempLat

namespace
{
  // 1D MemoryToolBox is rejected at compile time in MPI builds (memorytoolbox.h static_assert).
  TempLat::TDDContainer<TempLat::LazyFourierViewTester<double, 2>> test2;
  TempLat::TDDContainer<TempLat::LazyFourierViewTester<double, 3>> test3;

#ifdef HAVE_FFTFLOAT
  TempLat::TDDContainer<TempLat::LazyFourierViewTester<float, 2>> test2f;
  TempLat::TDDContainer<TempLat::LazyFourierViewTester<float, 3>> test3f;
#endif
} // namespace
