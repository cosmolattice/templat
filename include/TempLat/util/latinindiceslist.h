#ifndef TEMPLAT_UTIL_LATININDICESLIST_H
#define TEMPLAT_UTIL_LATININDICESLIST_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/util/rangeiteration/for_in_range.h"
#include "TempLat/util/rangeiteration/tagliteral.h"
#include "TempLat/lattice/field/assignablefieldcollection.h"

namespace TempLat
{

  /** @brief A class which an enumerated list which starts from 1 rather than 0.
   *
   *
   * Unit test: ctest -R test-latinindiceslist
   **/
  template <class... Args> class LatinIndicesList : public AssignableCollectionBase<LatinIndicesList<Args...>, Args...>
  {
  public:
    // Put public methods here. These should change very little over time.
    using AssignableCollectionBase<LatinIndicesList<Args...>, Args...>::operator=;

    LatinIndicesList(Args... args) : AssignableCollectionBase<LatinIndicesList<Args...>, Args...>(args...) {}

    template <int N> auto getComp(Tag<N> t) { return (*this)(t + 1_c); }

    static constexpr int SHIFTIND = 1;
    using Getter = GetComponent;

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
  };

  template <class... Args> auto make_latinindices_list(Args... args) { return LatinIndicesList<Args...>(args...); }

} // namespace TempLat

#endif
