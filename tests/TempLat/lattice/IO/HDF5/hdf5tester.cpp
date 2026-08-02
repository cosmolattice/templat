
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020
#ifdef HAVE_HDF5
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/IO/HDF5/filesaverhdf5.h"
#include "TempLat/lattice/IO/HDF5/fileloaderhdf5.h"
#include "TempLat/lattice/algebra/helpers/getvectorcomponent.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/util/almostequal.h"

namespace TempLat
{

  struct HDF5Tester {
    static void Test(TDDAssertion &tdd);
  };

  namespace
  {
    /** @brief Read a whole dataset back with plain HDF5, bypassing FileLoaderHDF5.
     *
     *  The loader has no sparse support at all -- it assumes the dataset covers the whole lattice
     *  -- so it cannot check a setLimits save. Reading the file directly also means the saver is
     *  being checked against the file format rather than against its own inverse.
     */
    std::vector<double> readDatasetDirect(const std::string &file, const std::string &name,
                                          std::vector<hsize_t> &dims)
    {
      hid_t f = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      hid_t d = H5Dopen2(f, name.c_str(), H5P_DEFAULT);
      hid_t s = H5Dget_space(d);

      dims.resize(H5Sget_simple_extent_ndims(s));
      H5Sget_simple_extent_dims(s, dims.data(), nullptr);

      size_t total = 1;
      for (auto n : dims)
        total *= n;

      std::vector<double> data(total);
      H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

      H5Sclose(s);
      H5Dclose(d);
      H5Fclose(f);
      return data;
    }

    /** @brief Every rank reads the file, so wait until every rank has finished writing it. */
    void barrier()
    {
#ifdef HAVE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
  } // namespace

  void HDF5Tester::Test(TDDAssertion &tdd)
  {
    // Existing test: save/load field round-trip
    {
      FileSaverHDF5 fs;
      FileLoaderHDF5 fl;

      const device::Idx nGrid = 16, nGhost = 1;
      auto toolBox = MemoryToolBox<3>::makeShared(nGrid, nGhost);

      Field<double, 3> phi("phi", toolBox);
      SpatialCoordinate<3> coords(toolBox);
      auto x = coords(1_c);
      auto y = coords(2_c);
      auto z = coords(3_c);
      auto local_idx = x * nGrid * nGrid + y * nGrid + z;
      phi = local_idx + 42.0;

      fs.create("./FILE.h5");
      fs.save(phi);
      fs.save(x * 1.);
      fs.save(0.45, "aDot");
      fs.close();

      Field<double, 3> psi("phi", toolBox);
      double aDot = 0;

      fl.open("./FILE.h5");
      fl.load(psi);
      fl.load(aDot, "aDot");
      fl.close();

      tdd.verify(AlmostEqual(aDot, 0.45));

      {
        auto localView = psi.getLocalNDHostView();
        bool all_correct = true;
        for (size_t i = 0; i < localView.extent(0); ++i)
          for (size_t j = 0; j < localView.extent(1); ++j)
            for (size_t k = 0; k < localView.extent(2); ++k) {
              all_correct &=
                  (AlmostEqual(localView(i, j, k), 42.0 + local_idx.eval(i + nGhost, j + nGhost, k + nGhost)));
              if (!AlmostEqual(localView(i, j, k), 42.0 + local_idx.eval(i + nGhost, j + nGhost, k + nGhost))) {
                std::cout << "Error at " << i << " " << j << " " << k << " got " << localView(i, j, k) << " expected "
                          << 42.0 + local_idx.eval(i + nGhost, j + nGhost, k + nGhost) << std::endl;
              }
            }
        tdd.verify(all_correct);
      }

      // The expression save above (fs.save(x * 1.)) used to be written and never looked at. It is
      // the only path left after saveBlock dropped the field fast path, so check its content: the
      // configuration-space spatial coordinate of index c is c itself (the sign conversion
      // midpoint is nGrid, not nGrid/2), so the dataset must be the dimension-0 index.
      barrier();
      {
        std::vector<hsize_t> dims;
        const auto raw = readDatasetDirect("./FILE.h5", "/x_0 * 1", dims);
        bool ok = (dims.size() == 3 && dims[0] == (hsize_t)nGrid && dims[1] == (hsize_t)nGrid &&
                   dims[2] == (hsize_t)nGrid);
        tdd.verify(ok, "Expression dataset has the full lattice extent");
        for (device::Idx i = 0; ok && i < nGrid; ++i)
          for (device::Idx j = 0; j < nGrid; ++j)
            for (device::Idx k = 0; k < nGrid; ++k)
              ok = ok && AlmostEqual(raw[(i * nGrid + j) * nGrid + k], (double)i);
        tdd.verify(ok, "Expression save writes the x coordinate at every lattice site");
      }
    }

    // Non-cubic lattice. Every other test here uses 16 in all three dimensions, which hides any
    // dimension swap or stride mix-up completely: with equal extents a transposed write still
    // lands inside the dataset and still round-trips through the equally transposed read.
    {
      FileSaverHDF5 fs;
      FileLoaderHDF5 fl;

      const device::Idx nGhost = 1;
      const device::IdxArray<3> nGrid{{8, 16, 32}};
      auto toolBox = MemoryToolBox<3>::makeShared(nGrid, nGhost);

      Field<double, 3> phi("phi", toolBox);
      SpatialCoordinate<3> coords(toolBox);
      auto x = coords(1_c);
      auto y = coords(2_c);
      auto z = coords(3_c);
      auto global_idx = x * (nGrid[1] * nGrid[2]) + y * nGrid[2] + z;
      phi = global_idx + 0.5;

      fs.create("./FILE_noncubic.h5");
      fs.save(phi);
      fs.close();

      Field<double, 3> psi("phi", toolBox);
      fl.open("./FILE_noncubic.h5");
      fl.load(psi);
      fl.close();

      {
        auto localView = psi.getLocalNDHostView();
        bool all_correct = true;
        for (size_t i = 0; i < localView.extent(0); ++i)
          for (size_t j = 0; j < localView.extent(1); ++j)
            for (size_t k = 0; k < localView.extent(2); ++k)
              all_correct &= AlmostEqual(localView(i, j, k),
                                         0.5 + global_idx.eval(i + nGhost, j + nGhost, k + nGhost));
        tdd.verify(all_correct, "Non-cubic 8x16x32 lattice round-trips through save and load");
      }

      barrier();
      {
        std::vector<hsize_t> dims;
        const auto raw = readDatasetDirect("./FILE_noncubic.h5", "/phi(x)", dims);
        bool ok = (dims.size() == 3 && dims[0] == 8 && dims[1] == 16 && dims[2] == 32);
        tdd.verify(ok, "Non-cubic dataset keeps the lattice dimension order 8x16x32");
        for (device::Idx i = 0; ok && i < 8; ++i)
          for (device::Idx j = 0; j < 16; ++j)
            for (device::Idx k = 0; k < 32; ++k)
              ok = ok && AlmostEqual(raw[(i * 16 + j) * 32 + k], 0.5 + (double)((i * 16 + j) * 32 + k));
        tdd.verify(ok, "Non-cubic dataset holds the right value at every index");
      }
    }

    // Sparse (subvolume) save. There was no test of setLimits at all, even though it changes the
    // saved geometry in every dimension at once.
    {
      const device::Idx nGrid = 16, nGhost = 1;
      auto toolBox = MemoryToolBox<3>::makeShared(nGrid, nGhost);

      Field<double, 3> phi("phi", toolBox);
      SpatialCoordinate<3> coords(toolBox);
      auto x = coords(1_c);
      auto y = coords(2_c);
      auto z = coords(3_c);
      auto global_idx = x * nGrid * nGrid + y * nGrid + z;
      phi = global_idx * 1.0;

      // Step divides (up - down) in every dimension: extents 6, 16, 7.
      {
        FileSaverHDF5 fs;
        const std::vector<device::Idx> down{2, 0, 1}, up{14, 16, 15}, step{2, 1, 2};
        fs.setLimits(down, up, step, 3);
        fs.create("./FILE_sparse.h5");
        fs.save(phi);
        fs.close();

        barrier();
        std::vector<hsize_t> dims;
        const auto raw = readDatasetDirect("./FILE_sparse.h5", "/phi(x)", dims);
        bool ok = (dims.size() == 3 && dims[0] == 6 && dims[1] == 16 && dims[2] == 7);
        tdd.verify(ok, "Sparse dataset extent is (up-down)/step in every dimension");
        for (device::Idx i = 0; ok && i < 6; ++i)
          for (device::Idx j = 0; j < 16; ++j)
            for (device::Idx k = 0; k < 7; ++k) {
              const double expected = (double)((2 + 2 * i) * nGrid * nGrid + j * nGrid + (1 + 2 * k));
              ok = ok && AlmostEqual(raw[(i * 16 + j) * 7 + k], expected);
            }
        tdd.verify(ok, "Sparse save samples down + step*index in every dimension");
      }

      // Step does NOT divide (up - down) in dimension 0: 15 - 2 = 13 is not a multiple of 2, so
      // the extent floors to 6 and coordinate 14 -- which would need index 6 -- is dropped. The
      // old recursion dropped it too, by handing H5Sselect_hyperslab an out-of-range offset and
      // ignoring the failure; here it is dropped deliberately and silently.
      {
        FileSaverHDF5 fs;
        const std::vector<device::Idx> down{2, 0, 1}, up{15, 16, 15}, step{2, 1, 2};
        fs.setLimits(down, up, step, 3);
        fs.create("./FILE_sparse_odd.h5");
        fs.save(phi);
        fs.close();

        barrier();
        std::vector<hsize_t> dims;
        const auto raw = readDatasetDirect("./FILE_sparse_odd.h5", "/phi(x)", dims);
        bool ok = (dims.size() == 3 && dims[0] == 6 && dims[1] == 16 && dims[2] == 7);
        tdd.verify(ok, "Non-divisible sparse extent floors to 6");
        for (device::Idx i = 0; ok && i < 6; ++i)
          for (device::Idx j = 0; j < 16; ++j)
            for (device::Idx k = 0; k < 7; ++k) {
              const double expected = (double)((2 + 2 * i) * nGrid * nGrid + j * nGrid + (1 + 2 * k));
              ok = ok && AlmostEqual(raw[(i * 16 + j) * 7 + k], expected);
            }
        tdd.verify(ok, "Non-divisible sparse save holds the first floor((up-down)/step) samples");
      }

    }

    // NDim 1 and 2. Besides being the only coverage of a non-3D lattice, NDim 1 is the only case
    // that reaches the RangePolicy branch of getLocalKokkosPolicy rather than an MDRangePolicy.
    // It is therefore only reachable in a serial build: memorytoolbox.h static_asserts NDim >= 2
    // under MPI, because FFTW's MPI interface returns a NULL plan for a rank-1 R2C transform.
    {
      const device::Idx nGhost = 1;
#ifndef HAVE_MPI
      {
        const device::Idx nGrid = 16;
        auto toolBox = MemoryToolBox<1>::makeShared(nGrid, nGhost);
        Field<double, 1> phi("phi", toolBox);
        SpatialCoordinate<1> coords(toolBox);
        auto x = coords(1_c);
        phi = x * 2.0;

        FileSaverHDF5 fs;
        fs.create("./FILE_1d.h5");
        fs.save(phi);
        fs.close();

        barrier();
        std::vector<hsize_t> dims;
        const auto raw = readDatasetDirect("./FILE_1d.h5", "/phi(x)", dims);
        bool ok = (dims.size() == 1 && dims[0] == (hsize_t)nGrid);
        tdd.verify(ok, "1D dataset has the lattice extent");
        for (device::Idx i = 0; ok && i < nGrid; ++i)
          ok = ok && AlmostEqual(raw[i], 2.0 * (double)i);
        tdd.verify(ok, "1D save writes the right value at every site");
      }
#endif
      {
        const device::IdxArray<2> nGrid{{8, 16}};
        auto toolBox = MemoryToolBox<2>::makeShared(nGrid, nGhost);
        Field<double, 2> phi("phi", toolBox);
        SpatialCoordinate<2> coords(toolBox);
        auto x = coords(1_c);
        auto y = coords(2_c);
        phi = x * nGrid[1] + y;

        FileSaverHDF5 fs;
        fs.create("./FILE_2d.h5");
        fs.save(phi);
        fs.close();

        barrier();
        std::vector<hsize_t> dims;
        const auto raw = readDatasetDirect("./FILE_2d.h5", "/phi(x)", dims);
        bool ok = (dims.size() == 2 && dims[0] == 8 && dims[1] == 16);
        tdd.verify(ok, "2D dataset keeps the lattice dimension order 8x16");
        for (device::Idx i = 0; ok && i < 8; ++i)
          for (device::Idx j = 0; j < 16; ++j)
            ok = ok && AlmostEqual(raw[i * 16 + j], (double)(i * 16 + j));
        tdd.verify(ok, "2D save writes the right value at every site");
      }
    }

    // Time series test
    {
      HDF5File file;
      file.create("test_extendible.h5");
      bool amIRoot = true;
#ifdef HAVE_MPI
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      amIRoot = (rank == 0);
#endif

      auto group = file.createOrOpenGroup("av");
      auto dataset = group.createTimeSeries<int>("t", {0}, {10});

      for (int i = 0; i < 10; ++i) {
        dataset.push(i);
      }

      dataset.extend(10);
      dataset.flush(amIRoot);

      auto sizes = dataset.getSizes();
      tdd.verify(sizes[0] == 10, "Time series has correct size after flush");

      dataset.close();

      auto dataset2 = group.createTimeSeries<double>("s", {0, 4}, {10, 1});

      for (int i = 0; i < 10; ++i) {
        dataset2.push(std::vector<double>{i / 6.0, 8, 9.0 / (i + 1), 3});
      }
      dataset2.extend(10);
      dataset2.flush(amIRoot);

      dataset2.close();
      group.close();
      file.close();

      tdd.verify(true);
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::HDF5Tester> test;
}

#endif
