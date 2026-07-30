#ifndef TEMPLAT_UTIL_TDD_TDDCONTAINERBASE_H
#define TEMPLAT_UTIL_TDD_TDDCONTAINERBASE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019

#include "TempLat/util/tdd/tddmacros.h"

namespace TempLat
{
  /** @brief A class which provides the interface for TDDContainers.
   *
   * Unit test: ctest -R test-tddcontainerbase
   */
  class TDDContainerBase
  {
  public:
    // Put public methods here. These should change very little over time.
    virtual int Test() = 0;

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */

  public:
    template <typename TestObjectUnknownHere> static inline void Test(TestObjectUnknownHere &tdd);
  };
} // namespace TempLat

#endif
