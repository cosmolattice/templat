#ifndef TEMPLAT_PARALLEL_DEVICE_GUARD_H
#define TEMPLAT_PARALLEL_DEVICE_GUARD_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2025

#include "TempLat/parallel/device.h"

#ifdef DEVICE_KOKKOS

#include "TempLat/parallel/devices/kokkos/session/kokkos_guard.h"

#elif defined(DEVICE_STD)

#include "TempLat/parallel/devices/std/std_guard.h"

#else

static_assert(false, "Unknown device type.");

#endif

namespace TempLat
{
  using export_device_namespace::DeviceGuard;
} // namespace TempLat

#endif
