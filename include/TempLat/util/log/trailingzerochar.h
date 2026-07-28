#ifndef TEMPLAT_UTIL_LOG_TRAILINGZEROCHAR_H
#define TEMPLAT_UTIL_LOG_TRAILINGZEROCHAR_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019

namespace TempLat
{

  /** @brief Performed by the compiler, if the const char* is a string at compile time. */
  constexpr const char *TrailingZeroChar(const char *name) { return *name == '\0' ? name : TrailingZeroChar(name + 1); }

} // namespace TempLat

#endif
