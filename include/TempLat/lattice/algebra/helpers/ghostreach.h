#ifndef TEMPLAT_LATTICE_ALGEBRA_HELPERS_GHOSTREACH_H
#define TEMPLAT_LATTICE_ALGEBRA_HELPERS_GHOSTREACH_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>

#include "TempLat/util/exception.h"

namespace TempLat
{
  MakeException(GhostDepthTooSmall);

  /** @brief The offsets at which an expression reads a field, relative to the site being evaluated: along each
   *  axis, the smallest (lo) and largest (hi) offset.
   *
   *  The ghost walk (ConfirmGhosts) carries it from the root of an expression down to the fields. Non-local nodes
   *  add their stencil (a shift by s adds s to both ends, a central difference widens by [-1, +1]), every other
   *  node passes it on unchanged. A field reached with a range of depth d needs d ghost layers. Since the ghost
   *  exchange also fills corners, the axes are independent: shift<1>(shift<2>(phi)) has depth 1.
   *
   * Unit test: ctest -R test-ghostreach
   **/
  struct GhostReach {
    static constexpr size_t maxNDim = 8;

    std::array<int, maxNDim> lo{}, hi{};

    /** @brief Offset by s along axis. */
    constexpr GhostReach shifted(size_t axis, int s) const { return widened(axis, s, s); }

    /** @brief Extend by a stencil reading offsets [down, up] along axis. */
    constexpr GhostReach widened(size_t axis, int down, int up) const
    {
      GhostReach res = *this;
      res.lo[axis] += down;
      res.hi[axis] += up;
      return res;
    }

    /** @brief Extend by a stencil reading offsets [down, up] along each of the first nDim axes. */
    constexpr GhostReach widenedAll(size_t nDim, int down, int up) const
    {
      GhostReach res = *this;
      for (size_t i = 0; i < nDim; ++i)
        res = res.widened(i, down, up);
      return res;
    }

    /** @brief Smallest range containing both. */
    constexpr GhostReach &operator|=(const GhostReach &other)
    {
      for (size_t i = 0; i < maxNDim; ++i) {
        lo[i] = std::min(lo[i], other.lo[i]);
        hi[i] = std::max(hi[i], other.hi[i]);
      }
      return *this;
    }

    constexpr bool isZero() const { return depth() == 0; }

    /** @brief Number of ghost layers needed to read all offsets in the range. */
    constexpr int depth() const
    {
      int d = 0;
      for (size_t i = 0; i < maxNDim; ++i)
        d = std::max({d, hi[i], -lo[i]});
      return d;
    }

    /** @brief The ranges along the first nDim axes, e.g. "[-1, 1] x [0, 0] x [0, 2]". */
    std::string toString(size_t nDim) const
    {
      std::string res;
      for (size_t i = 0; i < nDim; ++i)
        res += (i ? " x [" : "[") + std::to_string(lo[i]) + ", " + std::to_string(hi[i]) + "]";
      return res;
    }

    friend constexpr bool operator==(const GhostReach &, const GhostReach &) = default;
  };
} // namespace TempLat

#endif
