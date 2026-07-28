#ifndef TEMPLAT_UTIL_CONDITIONALOUTPUT_CONDITIONALSAYSHORT_H
#define TEMPLAT_UTIL_CONDITIONALOUTPUT_CONDITIONALSAYSHORT_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019

#include "TempLat/util/log/log.h"

namespace TempLat
{
  /** @brief A class which only spits out to TTY if enabled == true.
   *
   *
   * Unit test: ctest -R test-conditionalsayshort
   */
  class ConditionalSayShort
  {
  public:
    // Put public methods here. These should change very little over time.
    ConditionalSayShort(bool enabled_) : enabled(enabled_) {}

    template <typename T> ConditionalSayShort &operator<<(const T &t)
    {
      if (enabled) {
        sayShort << t;
      }
      return *this;
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    bool enabled;
  };
} // namespace TempLat

#endif
