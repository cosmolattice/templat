#ifndef TEMPLAT_PARALLEL_DEVICE_ITERATION_H
#define TEMPLAT_PARALLEL_DEVICE_ITERATION_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2025

#include "TempLat/parallel/device.h"

#ifdef DEVICE_KOKKOS

#include "TempLat/parallel/devices/kokkos/kokkos_iteration.h"

#elif defined(DEVICE_STD)

#include "TempLat/parallel/devices/std/std_iteration.h"

#else

static_assert(false, "No device iteration backend selected.");

#endif

namespace TempLat::device::iteration
{
  using export_device_namespace::iteration::fence;
  using export_device_namespace::iteration::foreach;
  using export_device_namespace::iteration::reduce;

  using export_device_namespace::iteration::Max;
  using export_device_namespace::iteration::Min;
  using export_device_namespace::iteration::Prod;
  using export_device_namespace::iteration::Sum;
} // namespace TempLat::device::iteration

#endif
