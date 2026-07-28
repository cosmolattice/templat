/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
/** \file The main function for the testing of a single unit testing. Link against all object files that you want to
 * test. */

#include "TempLat/util/tdd/tdd.h"
#include "TempLat/session/sessionguard.h"

int main(int argc, char *argv[])
{
  TempLat::SessionGuard guard(argc, argv);
  return TempLat::TDDRegister::run();
};
