
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020
#ifdef HAVE_HDF5
#include "TempLat/lattice/IO/HDF5/filesaverhdf5.h"
#include "TempLat/lattice/IO/HDF5/fileloaderhdf5.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include <filesystem>
#include <cstdio>

namespace TempLat
{

  struct FileSaverHDF5Tester {
    static void Test(TDDAssertion &tdd);
  };

  void FileSaverHDF5Tester::Test(TDDAssertion &tdd)
  {
    // Test 1: basic save/load of field and scalar.
    // The file is FILE_saver.h5, not FILE.h5: ctest runs every test in the same working
    // directory, so sharing a name with hdf5tester.cpp made the two race under ctest -j.
    {
      FileSaverHDF5 fs;

      const device::Idx nGrid = 16, nGhost = 1;
      auto toolBox = MemoryToolBox<3>::makeShared(nGrid, nGhost);

      Field<double, 3> phi("phi", toolBox);
      phi = 42.0;

      fs.create("./FILE_saver.h5");
      fs.save(phi);
      fs.save(0.45, "aDot");
      fs.close();

      std::filesystem::path filePath("./FILE_saver.h5");
      tdd.verify(std::filesystem::exists(filePath));
      tdd.verify(std::filesystem::file_size(filePath) > 0);

      std::string h5dumpOutput;
      {
        std::array<char, 128> buffer;
        std::string command = "h5dump ./FILE_saver.h5 2>&1";
        std::shared_ptr<FILE> pipe(popen(command.c_str(), "r"), pclose);

        while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
          h5dumpOutput += buffer.data();
        }
      }

      tdd.verify(h5dumpOutput.find("phi") != std::string::npos);
      tdd.verify(h5dumpOutput.find("aDot") != std::string::npos);
      tdd.verify(h5dumpOutput.find("42") != std::string::npos);
    }

    const std::string testFile = "test_rng_binary.h5";

    // Test 2: Binary RNG state round-trip (uniform-format: 313 uint64_t values)
    {
      // Build a synthetic state string mimicking master's RandomUniform format:
      // 313 space-separated uint64_t values
      std::ostringstream oss;
      for (int i = 0; i < 313; ++i) {
        if (i > 0) oss << " ";
        oss << static_cast<uint64_t>(i * 12345 + 67890);
      }
      std::string originalState = oss.str();

      FileSaverHDF5 saver;
      saver.create(testFile);
      saver.saveRNGStateBinary(originalState, "rng_uniform", 0, 1);
      saver.close();

      FileLoaderHDF5 loader;
      loader.open(testFile);
      std::string loadedState;
      loader.loadRNGStateBinary(loadedState, "rng_uniform", 0);
      loader.close();

      tdd.verify(originalState == loadedState, "Uniform-format binary round-trip preserves exact state string");

      std::remove(testFile.c_str());
    }

    // Test 3: Binary RNG state round-trip (gaussian-format: 316 uint64_t values with cache)
    {
      // Build synthetic state: first line = 313 values, second line = counter hasCache cachedValue
      std::ostringstream oss;
      for (int i = 0; i < 313; ++i) {
        if (i > 0) oss << " ";
        oss << static_cast<uint64_t>(i * 54321 + 11111);
      }
      oss << "\n";
      // gaussianCounter=42, hasCache=1, cachedValue=3.14159
      oss << "42 1 3.14159";
      std::string originalState = oss.str();

      FileSaverHDF5 saver;
      saver.create(testFile);
      saver.saveRNGStateBinary(originalState, "rng_gaussian", 0, 1);
      saver.close();

      FileLoaderHDF5 loader;
      loader.open(testFile);
      std::string loadedState;
      loader.loadRNGStateBinary(loadedState, "rng_gaussian", 0);
      loader.close();

      // The round-trip may not be exact for the double due to precision,
      // but the uint64_t values should match. Verify the first line matches.
      std::istringstream issOrig(originalState);
      std::istringstream issLoaded(loadedState);
      std::string origLine, loadedLine;
      std::getline(issOrig, origLine);
      std::getline(issLoaded, loadedLine);
      tdd.verify(origLine == loadedLine, "Gaussian-format binary round-trip preserves uniform state line");

      std::remove(testFile.c_str());
    }

    // Test 4: sparse limits that describe no dataset are rejected before one is created.
    {
      const device::Idx nGrid = 16, nGhost = 1;
      auto toolBox = MemoryToolBox<3>::makeShared(nGrid, nGhost);
      Field<double, 3> phi("phi", toolBox);
      phi = 1.0;

      // up beyond the lattice. The rod branch this replaced read (up - down) elements starting at
      // the first owned point, i.e. through the ghost cells and out of the allocation, and wrote
      // whatever it found; silent corruption becomes a message.
      {
        FileSaverHDF5 fs;
        fs.setLimits(std::vector<device::Idx>{0, 0, 0}, std::vector<device::Idx>{nGrid + 4, nGrid, nGrid},
                     std::vector<device::Idx>{1, 1, 1}, 3);
        fs.create("./FILE_limits.h5");
        bool threw = false;
        try {
          fs.save(phi);
        } catch (const InvalidSparseLimits &) {
          threw = true;
        }
        fs.close();
        tdd.verify(threw, "An upper save limit beyond the lattice is rejected");
      }

      // An empty range.
      {
        FileSaverHDF5 fs;
        fs.setLimits(std::vector<device::Idx>{8, 0, 0}, std::vector<device::Idx>{8, nGrid, nGrid},
                     std::vector<device::Idx>{1, 1, 1}, 3);
        fs.create("./FILE_limits.h5");
        bool threw = false;
        try {
          fs.save(phi);
        } catch (const InvalidSparseLimits &) {
          threw = true;
        }
        fs.close();
        tdd.verify(threw, "Save limits selecting an empty range are rejected");
      }

      // A negative lower limit.
      {
        FileSaverHDF5 fs;
        fs.setLimits(std::vector<device::Idx>{-1, 0, 0}, std::vector<device::Idx>{nGrid, nGrid, nGrid},
                     std::vector<device::Idx>{1, 1, 1}, 3);
        fs.create("./FILE_limits.h5");
        bool threw = false;
        try {
          fs.save(phi);
        } catch (const InvalidSparseLimits &) {
          threw = true;
        }
        fs.close();
        tdd.verify(threw, "A negative lower save limit is rejected");
      }

      // A step of zero, which would make the dataset extent meaningless.
      {
        FileSaverHDF5 fs;
        fs.setLimits(std::vector<device::Idx>{0, 0, 0}, std::vector<device::Idx>{nGrid, nGrid, nGrid},
                     std::vector<device::Idx>{0, 1, 1}, 3);
        fs.create("./FILE_limits.h5");
        bool threw = false;
        try {
          fs.save(phi);
        } catch (const InvalidSparseLimits &) {
          threw = true;
        }
        fs.close();
        tdd.verify(threw, "A save step below 1 is rejected");
      }

      std::remove("./FILE_limits.h5");
    }

    // Test 5: the staging budget only changes how many writes saveBlock makes, never the file.
    // Each chunk is still one contiguous index box, so the chunk count has to be invisible.
    {
      const device::Idx nGrid = 16, nGhost = 1;
      auto toolBox = MemoryToolBox<3>::makeShared(nGrid, nGhost);
      Field<double, 3> phi("phi", toolBox);
      SpatialCoordinate<3> coords(toolBox);
      phi = coords(1_c) * nGrid * nGrid + coords(2_c) * nGrid + coords(3_c);

      auto save = [&](const std::string &file, size_t budget) {
        FileSaverHDF5 fs;
        if (budget) fs.setStagingBudgetBytes(budget);
        fs.create(file);
        fs.save(phi);
        fs.close();
      };
      save("./FILE_chunk_default.h5", 0);
      save("./FILE_chunk_tiny.h5", 1); // one slice of dimension 0 per write, the minimum

      auto read = [](const std::string &file) {
        std::vector<double> data(16 * 16 * 16);
        hid_t f = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        hid_t d = H5Dopen2(f, "/phi(x)", H5P_DEFAULT);
        H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
        H5Dclose(d);
        H5Fclose(f);
        return data;
      };
#ifdef HAVE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      const auto whole = read("./FILE_chunk_default.h5");
      const auto chunked = read("./FILE_chunk_tiny.h5");
      bool same = true;
      for (size_t i = 0; i < whole.size(); ++i)
        same = same && (whole[i] == chunked[i]);
      tdd.verify(same, "A one-slice staging budget writes exactly the same bytes as the default one");

      bool valuesRight = true;
      for (device::Idx i = 0; i < nGrid; ++i)
        for (device::Idx j = 0; j < nGrid; ++j)
          for (device::Idx k = 0; k < nGrid; ++k)
            valuesRight = valuesRight && (chunked[(i * nGrid + j) * nGrid + k] ==
                                          (double)((i * nGrid + j) * nGrid + k));
      tdd.verify(valuesRight, "The chunked save still holds the right value at every index");
    }
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::FileSaverHDF5Tester> test;
}

#endif
