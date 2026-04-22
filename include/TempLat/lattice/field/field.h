#ifndef TEMPLAT_LATTICE_FIELD_FIELD_H
#define TEMPLAT_LATTICE_FIELD_FIELD_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
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
   **/
  template <typename T, size_t _NDim = 0> class Field : private FieldNDimCheck<_NDim>, public ConfigView<T, _NDim>
  {
  public:
    static constexpr size_t NDim = _NDim;
    using value_type = T;

    using ConfigView<T, NDim>::mManager;

    Field(std::string name, device::memory::host_ptr<MemoryToolBox<NDim>> toolBox,
          LatticeParameters<T> pLatPar = LatticeParameters<T>())
        : ConfigView<T, NDim>(name, toolBox, pLatPar), mFourierView(*this)
    {
    }

    template <typename R> void operator=(R &&g) { ConfigView<T, NDim>::operator=(g); }

    void operator=(const Field<T, NDim> &other) { operator=(OneType() * other); }

    FourierView<T, NDim> &inFourierSpace()
    {
      mManager->confirmFourierSpace();
      return mFourierView;
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
    FourierView<T, NDim> mFourierView;
  };
} // namespace TempLat

#endif
