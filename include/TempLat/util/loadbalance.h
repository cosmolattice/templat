#ifndef TEMPLAT_UTIL_LOADBALANCE_H
#define TEMPLAT_UTIL_LOADBALANCE_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include <cinttypes>
#include <cstddef>
#include <vector>

namespace TempLat
{

  /** @brief A class which
   * given and a number of threads nth, return nth integer which reprensents the most even splitting of the number
   *
   *
   * Unit test: ctest -R test-loadbalance
   **/
  class LoadBalance
  {
  public:
    // Put public methods here. These should change very little over time.
    LoadBalance() {}

    inline static std::vector<ptrdiff_t> getConf(ptrdiff_t nPoints, ptrdiff_t nThreads)
    {
      if (nThreads <= 0) return {}; // no threads: nothing to distribute (and avoid a divide-by-zero)
      auto quotRem = std::imaxdiv(nPoints, nThreads);
      std::vector<ptrdiff_t> res;
      for (ptrdiff_t i = 0; i < nThreads; ++i) {
        res.push_back(i < quotRem.rem ? quotRem.quot + 1 : quotRem.quot);
      }
      return res;
    }
  };

} // namespace TempLat

#endif
