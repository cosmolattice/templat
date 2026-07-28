#ifndef TEMPLAT_UTIL_STRINGIFY_H
#define TEMPLAT_UTIL_STRINGIFY_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020

/** @brief A macro used to parse preprocessor variable as strings, credits https://stackoverflow.com/a/34252436.
 *
 *
 *
 * Unit test: ctest -R test-stringify
 **/
#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)

#endif
