#ifndef TEMPLAT_LATTICE_IO_HDF5_FILELOADERHDF5_H
#define TEMPLAT_LATTICE_IO_HDF5_FILELOADERHDF5_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2020

#ifdef HAVE_HDF5

#include <algorithm>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <vector>
#include "TempLat/lattice/IO/HDF5/helpers/blockgeometry.h"
#include "TempLat/lattice/IO/HDF5/helpers/hdf5file.h"
#include "TempLat/lattice/algebra/helpers/getstring.h"
#include "TempLat/lattice/algebra/helpers/getgetreturntype.h"
#include "TempLat/parameters/parameterparser.h"

#include "TempLat/parallel/device.h"
#include "TempLat/parallel/device_iteration.h"
#include "TempLat/parallel/device_memory.h"

namespace TempLat
{

  /** @brief A class which implements loading in pure HDF5.
   *
   *
   * Unit test: ctest -R test-fileloaderhdf5
   **/
  class FileLoaderHDF5
  {
  public:
    // Put public methods here. These should change very little over time.

    FileLoaderHDF5() = default;

    /** @brief Open a file for reading.
     *
     *  Defaults to ReadOnly, and deliberately so: a *loader* has no write API at all -- every
     *  method below is an H5Dread -- so opening H5F_ACC_RDWR only ever bought the risk of one.
     *  HDF5 marks the superblock while a file is open for write and clears it on a clean close,
     *  so a job killed mid-load (walltime, node failure, scancel) could leave a production
     *  checkpoint flagged "already open for write" and unopenable until `h5clear -s`, on the one
     *  artifact a run cannot regenerate. Read-only opens never set that mark, never move the
     *  file's mtime, and work against a checkpoint directory the user has chmod'ed read-only.
     *
     *  The flag stays settable for the rare caller that really does want to reopen a file it
     *  intends to modify through the raw handle.
     */
    void open(std::string fn, FileMode flag = ReadOnly) { mFile.open(fn, flag); }

#ifdef HAVE_MPI
    /** @brief Point collective I/O at the communicator the lattice is decomposed over.
     *  Must be called before open. See HDF5File::setComm. */
    void setComm(MPI_Comm comm) { mFile.setComm(comm); }
#endif
    void close() { mFile.close(); }
    void reset() { this->close(); }

    void load(ParameterParser &par)
    {
      mDataset = mFile.openDataset("Parameters");

      std::vector<std::string> parStr;
      char tmp[HDF5TypeConstant::FixedSizeStringLength];

      // get number of parameters stored.
      hid_t dspace = H5Dget_space(mDataset);
      const int ndims = H5Sget_simple_extent_ndims(dspace);
      std::vector<hsize_t> dims(ndims);
      H5Sget_simple_extent_dims(dspace, dims.data(), NULL);
      auto nElements = dims[0];

      for (size_t i = 0; i < nElements; ++i) {
        mDataset.readElement(tmp, std::vector<hsize_t>(1, i));
        parStr.emplace_back(tmp);
      }
      H5Sclose(dspace);
      par.addFromVector(parStr);
      mDataset.close();
    }

    template <typename R> void load(R &t, std::string name)
    { // used to store a number. The name is the one of the dataset which contains this number.
      mDataset = mFile.openDataset(name);
      mDataset.readElement(&t, std::vector<hsize_t>(1, 0));
      mDataset.close();
    }

    /**
     * @brief Load a double scalar from a named dataset
     * @param value The value to load into
     * @param name Dataset name
     *
     * This overload explicitly handles double
     */
    void load(double &value, const std::string &name)
    {
      // Use "/" prefix for root group (matching openDataset pattern)
      std::string fullName = "/" + name;
      auto dataset = H5Dopen2(mFile.getHandle(), fullName.c_str(), H5P_DEFAULT);

      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);

      H5Dclose(dataset);
    }

    /**
     * @brief Load a uint64_t scalar from a named dataset
     * @param value The value to load into
     * @param name Dataset name
     */
    void loadScalarU64(uint64_t &value, const std::string &name)
    {
      std::string fullName = "/" + name;
      auto dataset = H5Dopen2(mFile.getHandle(), fullName.c_str(), H5P_DEFAULT);

      H5Dread(dataset, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);

      H5Dclose(dataset);
    }

    /**
     * @brief Check whether a named dataset exists in the file
     * @param name Dataset name
     *
     * Lets callers branch on checkpoint layout without provoking an HDF5 error.
     */
    bool datasetExists(const std::string &name)
    {
      std::string fullName = "/" + name;
      return H5Lexists(mFile.getHandle(), fullName.c_str(), H5P_DEFAULT) > 0;
    }

    /**
     * @brief Get the dimensions of a named dataset
     * @param name Dataset name
     * @return One entry per rank of the dataspace
     *
     * Used to distinguish checkpoint RNG layouts: legacy per-rank datasets are
     * [nRanks, stateLen], where stateLen == 1 for the counter-based generator
     * and 313/316 for the pre-Philox Mersenne-Twister state.
     */
    std::vector<hsize_t> getDatasetDims(const std::string &name)
    {
      std::string fullName = "/" + name;
      auto dataset = H5Dopen2(mFile.getHandle(), fullName.c_str(), H5P_DEFAULT);
      auto filespace = H5Dget_space(dataset);

      const int ndims = H5Sget_simple_extent_ndims(filespace);
      std::vector<hsize_t> dims(ndims > 0 ? ndims : 0);
      if (ndims > 0) H5Sget_simple_extent_dims(filespace, dims.data(), nullptr);

      H5Sclose(filespace);
      H5Dclose(dataset);
      return dims;
    }

    /**
     * @brief Load every element of a 1-D double dataset
     * @param values Output: resized to the dataset extent
     * @param name Dataset name
     *
     * loadPerRank reads only this rank's element. This reads the whole array,
     * which is what a rank-count change needs in order to re-aggregate
     * per-rank acceptance counters written by a differently sized run.
     */
    void loadWholeArray(std::vector<double> &values, const std::string &name)
    {
      std::string fullName = "/" + name;
      auto dataset = H5Dopen2(mFile.getHandle(), fullName.c_str(), H5P_DEFAULT);
      auto filespace = H5Dget_space(dataset);

      hsize_t dims[1];
      H5Sget_simple_extent_dims(filespace, dims, nullptr);
      values.assign(static_cast<size_t>(dims[0]), 0.0);

      auto plist = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
      H5Pset_dxpl_mpio(plist, H5FD_MPIO_INDEPENDENT);
#endif
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist, values.data());

      H5Pclose(plist);
      H5Sclose(filespace);
      H5Dclose(dataset);
    }

    /**
     * @brief Load a string from a named dataset
     * @param str String to load into
     * @param name Dataset name
     *
     * Uses a large fixed-size buffer (8KB) to accommodate RNG states (~5KB)
     */
    void load(std::string &str, const std::string &name)
    {
      constexpr size_t LargeStringLength = 16384; // 16KB for combined RNG states

      std::vector<char> buffer(LargeStringLength, 0);

      // Use "/" prefix for root group (matching openDataset pattern)
      std::string fullName = "/" + name;
      auto dataset = H5Dopen2(mFile.getHandle(), fullName.c_str(), H5P_DEFAULT);
      auto dtype = H5Dget_type(dataset);
      size_t typeSize = H5Tget_size(dtype);

      // Create memory type matching the file type size
      auto memtype = H5Tcopy(H5T_C_S1);
      H5Tset_size(memtype, typeSize);

      H5Dread(dataset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());

      H5Tclose(memtype);
      H5Tclose(dtype);
      H5Dclose(dataset);

      str = std::string(buffer.data());
    }

    /**
     * @brief Load per-rank string data from a shared dataset (parallel HDF5 safe)
     * @param str String to load into
     * @param name Dataset name (shared across all ranks)
     * @param mpiRank This rank's index
     *
     * Reads from a dataset of size [nRanks] where each rank reads its element.
     */
    void loadPerRank(std::string &str, const std::string &name, int mpiRank)
    {
      constexpr size_t LargeStringLength = 16384;

      std::vector<char> buffer(LargeStringLength, 0);

      std::string fullName = "/" + name;
      auto dataset = H5Dopen2(mFile.getHandle(), fullName.c_str(), H5P_DEFAULT);
      auto dtype = H5Dget_type(dataset);
      size_t typeSize = H5Tget_size(dtype);

      // Create memory type matching file type size
      auto memtype = H5Tcopy(H5T_C_S1);
      H5Tset_size(memtype, typeSize);

      // Select hyperslab for this rank's element
      auto filespace = H5Dget_space(dataset);
      hsize_t start[1] = {static_cast<hsize_t>(mpiRank)};
      hsize_t count[1] = {1};
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, count, nullptr);

      // Memory space for single element
      hsize_t memDims[1] = {1};
      auto memspace = H5Screate_simple(1, memDims, nullptr);

      // Read with independent I/O
      auto plist = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
      H5Pset_dxpl_mpio(plist, H5FD_MPIO_INDEPENDENT);
#endif
      H5Dread(dataset, memtype, memspace, filespace, plist, buffer.data());

      H5Pclose(plist);
      H5Sclose(memspace);
      H5Sclose(filespace);
      H5Tclose(memtype);
      H5Tclose(dtype);
      H5Dclose(dataset);

      str = std::string(buffer.data());
    }

    /**
     * @brief Load per-rank double data from a shared dataset (parallel HDF5 safe)
     * @param value The value to load into
     * @param name Dataset name (shared across all ranks)
     * @param mpiRank This rank's index
     */
    void loadPerRank(double &value, const std::string &name, int mpiRank)
    {
      std::string fullName = "/" + name;
      auto dataset = H5Dopen2(mFile.getHandle(), fullName.c_str(), H5P_DEFAULT);

      // Select hyperslab for this rank's element
      auto filespace = H5Dget_space(dataset);
      hsize_t start[1] = {static_cast<hsize_t>(mpiRank)};
      hsize_t count[1] = {1};
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, count, nullptr);

      // Memory space for single element
      hsize_t memDims[1] = {1};
      auto memspace = H5Screate_simple(1, memDims, nullptr);

      // Read with independent I/O
      auto plist = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
      H5Pset_dxpl_mpio(plist, H5FD_MPIO_INDEPENDENT);
#endif
      H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, plist, &value);

      H5Pclose(plist);
      H5Sclose(memspace);
      H5Sclose(filespace);
      H5Dclose(dataset);
    }

    /**
     * @brief Load RNG state from binary uint64_t array and reconstruct text string
     * @param textState Output: reconstructed text state for loadState()
     * @param name Dataset name
     * @param mpiRank This rank's index
     */
    void loadRNGStateBinary(std::string &textState, const std::string &name, int mpiRank)
    {
      std::string fullName = "/" + name;
      auto dataset = H5Dopen2(mFile.getHandle(), fullName.c_str(), H5P_DEFAULT);
      auto filespace = H5Dget_space(dataset);

      hsize_t dims[2];
      H5Sget_simple_extent_dims(filespace, dims, nullptr);
      size_t stateSize = dims[1];

      hsize_t start[2] = {static_cast<hsize_t>(mpiRank), 0};
      hsize_t count[2] = {1, stateSize};
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, count, nullptr);

      hsize_t memDims[1] = {stateSize};
      auto memspace = H5Screate_simple(1, memDims, nullptr);
      std::vector<uint64_t> binaryState(stateSize);

      auto plist = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
      H5Pset_dxpl_mpio(plist, H5FD_MPIO_INDEPENDENT);
#endif
      H5Dread(dataset, H5T_NATIVE_UINT64, memspace, filespace, plist, binaryState.data());

      H5Pclose(plist);
      H5Sclose(memspace);
      H5Dclose(dataset);
      H5Sclose(filespace);

      std::ostringstream oss;

      if (stateSize == 313) {
        for (size_t i = 0; i < stateSize; ++i) {
          if (i > 0) oss << " ";
          oss << binaryState[i];
        }
      } else if (stateSize == 316) {
        for (size_t i = 0; i < 313; ++i) {
          if (i > 0) oss << " ";
          oss << binaryState[i];
        }
        oss << "\n";
        oss << binaryState[313] << " " << binaryState[314] << " ";
        if (binaryState[314]) {
          double cachedValue;
          std::memcpy(&cachedValue, &binaryState[315], sizeof(double));
          oss << std::setprecision(17) << cachedValue;
        }
      } else {
        for (size_t i = 0; i < stateSize; ++i) {
          if (i > 0) oss << " ";
          oss << binaryState[i];
        }
      }

      textState = oss.str();
    }

    template <typename R> void load(R r)
    {
      mDataset = mFile.openDataset(GetString::get(r));
      loadBlock(r);
      mDataset.close();
    }

    /** @brief Read the block of the dataset this rank owns, with ONE H5Dread per chunk.
     *
     *  The mirror of FileSaverHDF5::saveBlock, and the same reason: the recursion this replaces
     *  issued one H5Dread per "rod" along the last lattice axis, so a 64^3 restart with 19 fields
     *  cost 77824 independent MPI-IO reads for 40 MB.
     *
     *  The loader has no sparse support -- it assumes the dataset covers the lattice one to one --
     *  so the geometry is the degenerate case of the saver's: ini = 0, end = nGrid, step = 1.
     *
     *  Public for the same reason saveBlock is: nvcc's extended-lambda restrictions do not allow
     *  the DEVICE_LAMBDA below inside a private member function.
     */
    template <typename R> void loadBlock(R r)
    {
      auto toolBox = r.getToolBox();
      constexpr size_t NDim = std::decay_t<decltype(*toolBox)>::NDim;
      using vType = typename GetGetReturnType<R>::type;

      const auto geo = computeDenseBlockGeometry<NDim>(toolBox->mLayouts.getConfigSpaceStarts(), // Local mpi offset.
                                                       toolBox->mLayouts.getConfigSpaceSizes(),  // Local mpi sizes.
                                                       toolBox->mNGridPointsVec);
      if (geo.isEmpty()) return;

      const auto mLayout = toolBox->mLayouts.getConfigSpaceLayout();

      device::IdxArray<NDim> memFirst = geo.first;
      device::apply([&](auto... idx) { mLayout.putMemoryIndexFromSpatialLocationInto(memFirst, idx...); }, memFirst);

      // Row-major strides of the staging buffer, matching H5Screate_simple's C ordering.
      device::IdxArray<NDim> bstride{};
      bstride[NDim - 1] = 1;
      for (size_t d = NDim - 1; d > 0; --d)
        bstride[d - 1] = bstride[d] * geo.count[d];
      const device::Idx sliceElems = (NDim > 1) ? bstride[0] : 1;

      const device::Idx maxSlices =
          std::max<device::Idx>(1, (device::Idx)(mStagingBudgetBytes / std::max<size_t>(1, sliceElems * sizeof(vType))));

      const device::IdxArray<NDim> starts{}; // foreach's third argument is an EXTENT, not a stop.
      const device::IdxArray<NDim> stepMem = geo.step;
      auto view = r.getView();

      for (device::Idx c0 = 0; c0 < geo.count[0]; c0 += maxSlices) {
        const device::Idx n0 = std::min(maxSlices, geo.count[0] - c0);

        device::IdxArray<NDim> extents = geo.count;
        device::IdxArray<NDim> offsets = geo.offset;
        device::IdxArray<NDim> memBase = memFirst;
        extents[0] = n0;
        offsets[0] += c0;
        memBase[0] += c0 * geo.step[0];

        const device::Idx nElems = n0 * sliceElems;

        std::vector<vType> host(nElems);
        mDataset.readSlices(host, extents, offsets);

        // Rank 1 keeps the staging view contiguous, so copyHostToDevice takes its fast path.
        device::memory::NDView<vType, 1> device_buf("LoadBlockBuffer", nElems);
        device::memory::copyHostToDevice(host.data(), device_buf);

        auto functor = DEVICE_LAMBDA(device::IdxArray<NDim> idx)
        {
          device::IdxArray<NDim> mem{};
          device::Idx linear = 0;
          for (size_t d = 0; d < NDim; ++d) {
            mem[d] = memBase[d] + idx[d] * stepMem[d];
            linear += idx[d] * bstride[d];
          }
          device::apply([&](const auto &...i) { view(i...) = device_buf(linear); }, mem);
        };
        device::iteration::foreach<NDim>("LoadBlockBufferScattering", starts, extents, functor);
      }
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    HDF5File mFile;
    HDF5Dataset mDataset;
    /** @brief Cap on the buffer loadBlock stages a field through; see
     *  FileSaverHDF5::setStagingBudgetBytes. */
    size_t mStagingBudgetBytes = 256ull << 20;
  };
} // namespace TempLat

#endif

#endif
