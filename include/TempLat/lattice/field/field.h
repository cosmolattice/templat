#ifndef TEMPLAT_LATTICE_FIELD_FIELD_H
#define TEMPLAT_LATTICE_FIELD_FIELD_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Franz R. Sattler, Year: 2025

#include "TempLat/lattice/algebra/constants/halftype.h"
#include "TempLat/lattice/algebra/constants/onetype.h"
#include "TempLat/lattice/algebra/constants/zerotype.h"
#include "TempLat/lattice/field/views/fieldviewconfig.h"
#include "TempLat/lattice/field/views/fieldviewfourier.h"
#include "TempLat/parallel/device.h"

namespace TempLat
{
  template <size_t NDimCheck> struct FieldNDimCheck {
    static_assert(NDimCheck > 0, "NDim template parameter is required. Use e.g. Field<double, 3>.");
  };

  /** @brief A class which is a classical field on your n-dimensional equisized grid.
   * You use it as a scalar field, a vector component, whatever.
   * Template parameter is your type of floating point precision: float or double. Default: double.
   *
   *  Implements a get method, and is hence suitable for all algebra.
   *
   * Unit test: ctest -R test-field
   *
   * @vocab-summary A real scalar field on the lattice, and the type that owns the memory. Assigning to it
   * evaluates an expression across every site; passing it into one makes it a leaf.
   * @vocab-signature Field<T, NDim> phi("phi", toolBox);
   **/
  template <typename T, size_t _NDim = 0> class Field : private FieldNDimCheck<_NDim>, public ConfigView<T, _NDim>
  {
  public:
    static constexpr size_t NDim = _NDim;
    using value_type = T;

    using ConfigView<T, NDim>::mManager;

    Field(std::string name, device::memory::host_ptr<MemoryToolBox<NDim>> toolBox,
          LatticeParameters<T> pLatPar = LatticeParameters<T>())
        : ConfigView<T, NDim>(name, toolBox, pLatPar)
    {
    }

    template <typename R> void operator=(R &&g) { ConfigView<T, NDim>::operator=(g); }

    void operator=(const Field<T, NDim> &other) { operator=(OneType() * other); }

    /** @brief The Fourier-space view, materialised on first use.
     *
     * confirmFourierSpace() begins with allocate(), so materialising after it is safe, and
     * the view cannot go stale afterwards: MemoryManager::allocate latches mAllocated and is
     * the only caller of MemoryBlock::allocate, so the data pointer and the extents are fixed
     * for the manager's lifetime, and FourierView caches no host view -- its three host
     * accessors re-derive from mManager on every call.
     *
     * Still returns FourierView& and not the handle: `f.inFourierSpace() = expr` and
     * .setZeroMode(...) depend on that.
     *
     * The lazy write is thread-safe ONLY because there is no host-side threading: there are
     * zero host `#pragma omp` regions in templat, CosmoInterface or MCInterface, and Kokkos
     * parallel regions run device lambdas, which never call this. If host-side OpenMP is ever
     * put around field access, this becomes a genuine data race. Not `mutable`, and this must
     * stay non-const, for the same reason: a const overload would make the race reachable
     * through a shared const&.
     */
    FourierView<T, NDim> &inFourierSpace()
    {
      mManager->confirmFourierSpace(); // allocates first
      // Constructed here rather than handed to host_ptr's variadic constructor: that one is
      // constrained by `requires requires(ARGS... args) { T(args...); }`, evaluated in
      // host_ptr's access context, and host_ptr is not a friend of FourierView. The
      // constraint is unsatisfied, the variadic overload drops out and nothing else is
      // viable. Inside Field the private constructor IS reachable -- fieldviewfourier.h
      // declares every Field specialisation a friend -- so build it and hand the raw pointer
      // to the unconstrained host_ptr(T*).
      if (mFourierView.get() == nullptr) mFourierView = new FourierView<T, NDim>(*this);
      return *mFourierView;
    }

    template <typename S>
      requires(!std::is_same_v<Field<T, NDim>, S>)
    auto d(const S &other) const
    {
      return ZeroType();
    }
    /** @brief The real overlord: is it a Field, then we must compare. */
    device::Idx d(const Field<T, NDim> &other) const { return *this == other ? 1 : 0; }

    friend bool operator==(const Field<T, NDim> &a, const Field<T, NDim> &b) { return a.mManager == b.mManager; }

    template <typename S>
      requires std::is_same_v<Field<T, NDim>, S>
    auto d(const S &other) const
    {
      return OneType();
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */

    /** @brief Behind a handle, not by value.
     *
     * FourierView derives from AbstractField, so held by value it contributed both its own
     * bytes and a second host_ptr<MemoryManager> reference count to every copy of every
     * Field -- and Fields are copied by value into every node of every expression tree.
     * Configuration-space expressions, which is all of MC and cooling, never touch it.
     * A field that never enters Fourier space now pays literally nothing for it: host_ptr's
     * copy constructor skips the atomic increment when the source is null and its destructor
     * returns early.
     *
     * It must NOT move into the MemoryManager: FourierView holds a host_ptr to the manager,
     * so that would close a reference-count cycle and leak every field's memory. As a handle
     * held by Field the ownership graph stays a diamond -- the heap view dies with the last
     * Field copy holding it, releasing its manager handle at the same instant it does today.
     *
     * Sharing between Field copies is unobservable: FourierView writes its members only in
     * its private constructor, holds no FFT plans (those are in MemoryToolBox) and no
     * space-state flags (those are in MemoryManager), and its implicit copy assignment is
     * deleted by `const LatticeParameters<T> latPar`, so nobody can overwrite a shared view
     * through the returned reference either. Two views built from the same Field are
     * bit-identical and nothing anywhere tests pointer identity.
     *
     * Never dereference this in device code. host_ptr's copy constructor nulls out under
     * DEVICE_REGION, so a Field captured into a KOKKOS_CLASS_LAMBDA carries a null handle on
     * the device -- harmless, exactly as for mManager today, but only while nothing reads it
     * there.
     */
    device::memory::host_ptr<FourierView<T, NDim>> mFourierView;
  };
} // namespace TempLat

#endif
