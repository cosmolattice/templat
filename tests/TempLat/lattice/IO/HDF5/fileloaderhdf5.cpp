/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2026
#ifdef HAVE_HDF5
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/IO/HDF5/filesaverhdf5.h"
#include "TempLat/lattice/IO/HDF5/fileloaderhdf5.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/util/almostequal.h"

#include <chrono>
#include <cstdio>
#include <filesystem>

namespace TempLat
{

  struct FileLoaderHDF5Tester {
    static void Test(TDDAssertion &tdd);
  };

  namespace
  {
    void barrier()
    {
#ifdef HAVE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    bool amIRoot()
    {
#ifdef HAVE_MPI
      int rank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      return rank == 0;
#else
      return true;
#endif
    }
  } // namespace

  /** @brief The loader must never need write access to the file it reads.
   *
   *  This is the regression test for a cooling stream having moved a production checkpoint's
   *  mtime just by loading it: FileLoaderHDF5::open forwarded no FileMode, so HDF5File::open
   *  took its ReadWrite default and issued H5Fopen(..., H5F_ACC_RDWR). Nothing was written, but
   *  RDWR marks the superblock until a clean close, so a job killed mid-load could leave the one
   *  artifact a run cannot regenerate flagged "already open for write".
   *
   *  A read-only file is the assertion that makes this impossible to regress: with RDWR the open
   *  below fails outright, and it also proves a checkpoint tree can be chmod'ed read-only.
   */
  void FileLoaderHDF5Tester::Test(TDDAssertion &tdd)
  {
    const std::string fileName = "./FILE_loader_readonly.h5";
    const std::filesystem::path filePath(fileName);

    const device::Idx nGrid = 8, nGhost = 1;
    auto toolBox = MemoryToolBox<3>::makeShared(nGrid, nGhost);

    {
      Field<double, 3> phi("phi", toolBox);
      phi = 17.5;

      FileSaverHDF5 fs;
      fs.create(fileName);
      fs.save(phi);
      fs.save(0.25, "aDot");
      fs.close();
    }
    barrier();

    // Strip every write bit, and back-date the file, so that an RDWR open both fails and would
    // be visible in the timestamp if it somehow succeeded.
    const auto backdated = std::filesystem::file_time_type::clock::now() - std::chrono::hours(1);
    if (amIRoot()) {
      std::filesystem::permissions(filePath, std::filesystem::perms::owner_write |
                                                 std::filesystem::perms::group_write |
                                                 std::filesystem::perms::others_write,
                                   std::filesystem::perm_options::remove);
      std::filesystem::last_write_time(filePath, backdated);
    }
    barrier();

    const auto mtimeBefore = std::filesystem::last_write_time(filePath);

    Field<double, 3> psi("phi", toolBox);
    psi = 0.0;
    double aDot = 0;

    {
      FileLoaderHDF5 fl;
      fl.open(fileName); // default FileMode must be ReadOnly
      fl.load(psi);
      fl.load(aDot, "aDot");
      fl.close();
    }
    barrier();

    tdd.verify(AlmostEqual(aDot, 0.25), "A write-protected file still loads its scalars");

    {
      auto localView = psi.getLocalNDHostView();
      bool allCorrect = true;
      for (size_t i = 0; i < localView.extent(0); ++i)
        for (size_t j = 0; j < localView.extent(1); ++j)
          for (size_t k = 0; k < localView.extent(2); ++k)
            allCorrect = allCorrect && AlmostEqual(localView(i, j, k), 17.5);
      tdd.verify(allCorrect, "A write-protected file still loads its fields");
    }

    const auto mtimeAfter = std::filesystem::last_write_time(filePath);
    tdd.verify(mtimeAfter == mtimeBefore, "Loading a file does not move its modification time");

    barrier();
    if (amIRoot()) {
      std::filesystem::permissions(filePath, std::filesystem::perms::owner_write,
                                   std::filesystem::perm_options::add);
      std::remove(fileName.c_str());
    }
    barrier();
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::FileLoaderHDF5Tester> test;
}

#endif
