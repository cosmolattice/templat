#ifndef TEMPLAT_UTIL_STATICIF_H
#define TEMPLAT_UTIL_STATICIF_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

namespace TempLat
{

/**
 * @vocab-summary Compile-time branch inside an expression: only the taken side is instantiated.
 * @vocab-signature IfElse(condition, ifExpr, elseExpr)
 **/
#define IfElse(condition, ifExpr, elseExpr)                                                                            \
  [&]() {                                                                                                              \
    if constexpr (condition) {                                                                                         \
      return ifExpr;                                                                                                   \
    } else {                                                                                                           \
      return elseExpr;                                                                                                 \
    }                                                                                                                  \
  }()
/**
 * @vocab-summary Compile-time branch whose untaken side is ZeroType, so it prunes out of the tree entirely.
 * @vocab-signature If(condition, ifExpr)
 **/
#define If(condition, ifExpr)                                                                                          \
  [&]() {                                                                                                              \
    if constexpr (condition) {                                                                                         \
      return ifExpr;                                                                                                   \
    } else {                                                                                                           \
      return ZeroType();                                                                                               \
    }                                                                                                                  \
  }()

#define IfElseStatement(condition, ifExpr, elseExpr)                                                                   \
  [&]() {                                                                                                              \
    if constexpr (condition) {                                                                                         \
      ifExpr;                                                                                                          \
    } else {                                                                                                           \
      elseExpr;                                                                                                        \
    }                                                                                                                  \
  }()
#define IfStatement(condition, ifExpr)                                                                                 \
  if constexpr (condition) {                                                                                           \
    ifExpr;                                                                                                            \
  }
} // namespace TempLat

#endif
