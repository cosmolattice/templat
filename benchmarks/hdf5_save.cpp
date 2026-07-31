/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2026

/** \file The cost of writing a field snapshot.
 *
 * Reproduces a production MC checkpoint: 19 datasets of a 64^3 lattice of doubles, 40 MB in
 * total, written to one HDF5 file. That shape is the reason this benchmark exists -- FileSaverHDF5
 * used to issue one H5Dwrite per "rod" (one contiguous run along the last lattice axis), so the
 * same 40 MB cost 19 * 64 * 64 = 77824 independent MPI-IO writes and took 17.1 s on a single GPU
 * rank, 2.3 MB/s against the ~270 MB/s the filesystem gives. The wall clock tracked rods per rank
 * rather than bytes, so it got *faster* with more ranks. Run this at 1, 4 and 32 ranks: the time
 * should now be flat-to-decreasing in the rank count, and the writes-per-rank line should read one
 * per dataset.
 *
 * Build with -DTEMPLAT_BENCH=ON, then
 *
 *     mpiexec -n <ranks> ./bench-hdf5_save
 *
 * Overridable through the environment: TEMPLAT_BENCH_NGRID (default 64), TEMPLAT_BENCH_NDATASETS
 * (default 19), TEMPLAT_BENCH_REPS (default 3), TEMPLAT_BENCH_FILE (default ./hdf5_save_bench.h5).
 */

#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "TempLat/session/sessionguard.h"
#include "TempLat/util/timer.h"

#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/lattice/IO/HDF5/filesaverhdf5.h"

namespace
{
  size_t fromEnv(const char *name, size_t fallback)
  {
    const char *value = std::getenv(name);
    if (value == nullptr) return fallback;
    const long parsed = std::strtol(value, nullptr, 10);
    return parsed > 0 ? (size_t)parsed : fallback;
  }
} // namespace

int main(int argc, char **argv)
{
  using namespace TempLat;

#ifndef HAVE_HDF5
  (void)argc;
  (void)argv;
  std::cout << "bench-hdf5_save needs -DHDF5=ON; nothing to measure.\n";
  return 0;
#else
  SessionGuard guard(argc, argv, false);

  constexpr size_t NDim = 3;
  const device::Idx nGrid = (device::Idx)fromEnv("TEMPLAT_BENCH_NGRID", 64);
  const size_t nDatasets = fromEnv("TEMPLAT_BENCH_NDATASETS", 19);
  const size_t nReps = fromEnv("TEMPLAT_BENCH_REPS", 3);
  const char *fileEnv = std::getenv("TEMPLAT_BENCH_FILE");
  const std::string file = fileEnv ? fileEnv : "./hdf5_save_bench.h5";

  auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, 1);
  toolBox->unsetVerbose();

  const int nRanks = toolBox->getNProcesses();
  const bool iAmRoot = toolBox->getMPIRank() == 0;

  std::vector<Field<double, NDim>> fields;
  fields.reserve(nDatasets);
  SpatialCoordinate<NDim> coords(toolBox);
  for (size_t i = 0; i < nDatasets; ++i) {
    fields.emplace_back("bench_" + std::to_string(i), toolBox);
    fields.back() = coords(1_c) * nGrid * nGrid + coords(2_c) * nGrid + coords(3_c) + (double)i;
  }

  const double megabytes = (double)nDatasets * (double)nGrid * (double)nGrid * (double)nGrid * sizeof(double) / 1e6;

  // The local sizes are the same on every rank of a regular decomposition, so counting on one
  // rank is enough to say how many writes a rank now issues -- and how many it used to.
  const auto localSizes = toolBox->mLayouts.getConfigSpaceSizes();
  size_t rodsPerRank = 1;
  for (size_t d = 0; d + 1 < NDim; ++d)
    rodsPerRank *= (size_t)localSizes[d];

  if (iAmRoot) {
    std::cout << "lattice " << nGrid << "^" << NDim << ", " << nDatasets << " datasets, " << std::fixed
              << std::setprecision(1) << megabytes << " MB, " << nRanks << " rank(s)\n"
              << "writes per rank: " << nDatasets << " (one per dataset); "
              << "the per-rod saver issued " << nDatasets * rodsPerRank << "\n";
  }

  for (size_t rep = 0; rep < nReps; ++rep) {
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    const Timer timer;

    FileSaverHDF5 saver;
    saver.create(file);
    for (auto &field : fields)
      saver.save(field);
    saver.close();

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    const double seconds = (double)timer.nanoseconds() * 1e-9;

    if (iAmRoot) {
      std::cout << "rep " << rep << ": " << std::fixed << std::setprecision(3) << seconds << " s, "
                << std::setprecision(1) << megabytes / seconds << " MB/s\n";
    }
  }

  if (iAmRoot) std::remove(file.c_str());
  return 0;
#endif
}
