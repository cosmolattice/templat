#ifndef TEMPLAT_UTIL_CONTAINSSPACE_H
#define TEMPLAT_UTIL_CONTAINSSPACE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019

#include <cctype>
#include <string>

namespace TempLat
{

  /** @brief A class which checks if a string contains a space.
   *
   * Unit test: ctest -R test-containsspace
   **/
  class ContainsSpace
  {
  public:
    // Put public methods here. These should change very little over time.

    static inline bool test(const std::string &input)
    {
      bool result = false;
      for (auto &&it : input)
        result = result || isspace(it);
      return result;
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    ContainsSpace() = default;
  };
} // namespace TempLat

#endif
