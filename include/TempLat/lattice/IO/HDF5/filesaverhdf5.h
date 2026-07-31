#ifndef TEMPLAT_LATTICE_IO_HDF5_FILESAVERHDF5_H
#define TEMPLAT_LATTICE_IO_HDF5_FILESAVERHDF5_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Franz R. Sattler, Year: 2025

#ifdef HAVE_HDF5

#include <algorithm>
#include <cstring>
#include <sstream>
#include <vector>
#include "TempLat/util/prettytostring.h"
#include "TempLat/lattice/algebra/helpers/ghostshunter.h"
#include "TempLat/lattice/algebra/helpers/confirmspace.h"
#include "TempLat/lattice/algebra/spacestateinterface.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/lattice/algebra/helpers/getgetreturntype.h"
#include "TempLat/lattice/algebra/helpers/getfloattype.h"
#include "TempLat/lattice/algebra/helpers/getstring.h"
#include "TempLat/parameters/parameterparser.h"

#include "TempLat/parallel/device.h"
#include "TempLat/parallel/device_iteration.h"

#include "TempLat/lattice/IO/HDF5/helpers/hdf5file.h"
#include "TempLat/lattice/IO/HDF5/helpers/blockgeometry.h"

namespace TempLat
{
  MakeException(StringIsTooLong);
  MakeException(InvalidSparseLimits);

  /** @brief A class which implements saving in pure HDF5.
   *
   *
   * Unit test: ctest -R test-filesaverhdf5
   **/
  class FileSaverHDF5
  {
  public:
    /* Put public methods here. These should change very little over time. */
    FileSaverHDF5() : mdown(10, 0), mstep(10, 1) {}

    void open(std::string fn) { mFile.open(fn); }

    void create(std::string fn, FileMode flag = Overwrite) { mFile.create(fn, flag); }

#ifdef HAVE_MPI
    /** @brief Point collective I/O at the communicator the lattice is decomposed over.
     *  Must be called before open/create. See HDF5File::setComm. */
    void setComm(MPI_Comm comm) { mFile.setComm(comm); }
#endif

    void close() { mFile.close(); }
    void reset() { this->close(); }

    template <typename R> void setLimits(R down, R up, R step, device::Idx NDim)
    {
      mdown.resize(NDim);
      mup.resize(NDim);
      mstep.resize(NDim);
      sparsesizes.resize(NDim);
      sparsesave = true;
      for (device::Idx i = 0; i < NDim; i++) {
        mdown[i] = down[i];
        mup[i] = up[i];
        mstep[i] = step[i];
        // Validation lives in save(), which is the first point that knows the lattice size -- but
        // a step of zero has to be sidestepped here, since it would divide by zero on the way.
        sparsesizes[i] = (mstep[i] >= 1) ? (mup[i] - mdown[i]) / mstep[i] : 0;
      }
    }

    /** @brief Cap on the buffer saveBlock stages a field into, per rank, before writing it.
     *
     *  saveBlock chunks the slowest dimension so that neither the device buffer nor its host copy
     *  exceeds this. Lower it if checkpointing a large local subdomain runs the device out of
     *  memory; raise it to trade memory for fewer, larger writes. The default is 256 MiB, well
     *  above any one field of a production lattice (2 MiB per rank per dataset at 64^3 doubles).
     */
    void setStagingBudgetBytes(size_t bytes) { mStagingBudgetBytes = bytes; }

    void save(ParameterParser &r)
    { // Conceptually, may be better as attributes? But nightmare to save vector of strings, did nt manage to do it in a
      // finite amount of time.
      std::ostringstream oss;

      std::vector<std::string> parStr;
      std::string tmp;
      for (auto x : r.getParams()) {
        if (x.second != "inf") {
          parStr.emplace_back(x.first + "=" + x.second);
        }
      }

      // Create a flat buffer of fixed-size strings
      std::vector<char> stringData(parStr.size() * HDF5TypeConstant::FixedSizeStringLength, 0);
      for (size_t i = 0; i < parStr.size(); ++i) {
        if (parStr[i].size() > HDF5TypeConstant::FixedSizeStringLength)
          throw StringIsTooLong("Well, that's a bit embarassing. One of your parameters contains too many characters "
                                "(the total string should be smaller than " +
                                std::to_string(HDF5TypeConstant::FixedSizeStringLength) +
                                " char by default, for our hdf5). If you managed to make HDF5 with variable string "
                                "length, please let us know! If you just want to change the hardcoded number, look in "
                                "the file TempLat/lattice/IO/HDF5/helpers/hdf5type.h .");
        std::strncpy(&stringData[i * HDF5TypeConstant::FixedSizeStringLength], parStr[i].c_str(),
                     HDF5TypeConstant::FixedSizeStringLength - 1);
      }

      // Create dataset and write directly with H5Dwrite using proper string type
      mDataset = mFile.createDataset<const char *>("Parameters", std::vector<hsize_t>(1, parStr.size()));
      HDF5Type<const char *> strtype;
      H5Dwrite(mDataset, strtype.type, H5S_ALL, H5S_ALL, H5P_DEFAULT, stringData.data());
      strtype.close();
      mDataset.close();
    }

    void save_attr(ParameterParser &r)
    { // Conceptually, may be better as attributes? But nightmare to save vector of strings, did nt manage to do it in a
      // finite amount of time.
      std::ostringstream oss;

      std::vector<std::string> parStr;
      std::string tmp;
      mDataset = mFile.createDataset<const char *>("Parameters", std::vector<hsize_t>(1, parStr.size()));
      for (auto x : r.getParams()) {
        if (x.second != "inf") {
          // parStr.emplace_back(x.first + "=" + x.second);
          mDataset.addAtribute(x.first, x.second);
        }
      }
      /* for(size_t i = 0; i < parStr.size(); ++i){
           if(parStr[i].size() > HDF5TypeConstant::FixedSizeStringLength) throw StringIsTooLong("Well, that's a bit
       embarassing. One of your parameters contains too many characters (the total string should be smaller than
       "+std::to_string(HDF5TypeConstant::FixedSizeStringLength)+" char by default, for our hdf5). If you managed to
       make HDF5 with variable string length, please let us know! If you just want to change the hardcoded number, look
       in the file TempLat/lattice/IO/HDF5/helpers/hdf5type.h .");
           mDataset.writeElement(parStr[i].c_str(),std::vector<hsize_t>(1,i));
       }*/
      mDataset.close();
    }

    template <typename R> void save(R r)
    { // used to store an entity directly to a dataset, using it's own name.
      using vType = GetGetReturnType<R>::type;
      ConfirmSpace::apply(r, r.getToolBox()->mLayouts.getConfigSpaceLayout(), SpaceStateType::Configuration);
      GhostsHunter::apply(r);
      if (!sparsesave)
        sparsesizes.assign(r.getToolBox()->mNGridPointsVec.begin(), r.getToolBox()->mNGridPointsVec.end());
      checkSparseLimits(r.getToolBox()->mNGridPointsVec, std::decay_t<decltype(*r.getToolBox())>::NDim);
      mDataset = mFile.createDataset<vType>(GetString::get(r), sparsesizes);
      saveBlock(r);
      mDataset.close();
    }

    template <typename R, typename T> void save(T t, R r, std::string name)
    { // used to store an entity in a time series. The name is the one of the group, data set labelled by t.
      using vType = GetGetReturnType<R>::type;
      ConfirmSpace::apply(r, r.getToolBox()->mLayouts.getConfigSpaceLayout(), SpaceStateType::Configuration);
      GhostsHunter::apply(r);
      if (!sparsesave)
        sparsesizes.assign(r.getToolBox()->mNGridPointsVec.begin(), r.getToolBox()->mNGridPointsVec.end());
      checkSparseLimits(r.getToolBox()->mNGridPointsVec, std::decay_t<decltype(*r.getToolBox())>::NDim);
      mDataset = mFile.createOrOpenGroup(name).createDataset<vType>(PrettyToString::get(t, 10), sparsesizes);
      saveBlock(r);
      mDataset.close();
    }

    template <typename R> void save(R t, std::string name)
    { // used to store a number. The name is the one of the dataset which contains this number.
      using vType = GetGetReturnType<R>::type;
      mDataset = mFile.createDataset<vType>(name, std::vector<hsize_t>(1, 1));
      // writeElement takes the value by value and passes &data to H5Dwrite itself; passing &t here would
      // deduce T = R* and store the pointer's address bits instead of the scalar.
      mDataset.writeElement(t, std::vector<hsize_t>(1, 0));
      mDataset.close();
    }

    /**
     * @brief Save a double scalar to a named dataset
     * @param value The value to save
     * @param name Dataset name
     *
     * This overload explicitly handles double without GetGetReturnType
     */
    void save(double value, const std::string &name)
    {
      // Use "/" prefix for root group (matching openDataset pattern)
      std::string fullName = "/" + name;
      hsize_t dims[1] = {1};
      auto dataspace = H5Screate_simple(1, dims, nullptr);
      auto dataset = H5Dcreate2(mFile.getHandle(), fullName.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT);

      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);

      H5Dclose(dataset);
      H5Sclose(dataspace);
    }

    /**
     * @brief Save a string to a named dataset
     * @param str The string to save
     * @param name Dataset name
     *
     * Uses a large fixed-size buffer (8KB) to accommodate RNG states (~5KB)
     */
    void save(const std::string &str, const std::string &name)
    {
      constexpr size_t LargeStringLength = 16384; // 16KB for combined RNG states

      if (str.size() >= LargeStringLength) {
        throw StringIsTooLong("String too long for HDF5 dataset '" + name + "': " + std::to_string(str.size()) +
                              " chars (max " + std::to_string(LargeStringLength) + ")");
      }

      std::vector<char> buffer(LargeStringLength, 0);
      std::strncpy(buffer.data(), str.c_str(), LargeStringLength - 1);

      // Create custom string type with larger size
      auto memtype = H5Tcopy(H5T_C_S1);
      H5Tset_size(memtype, LargeStringLength);

      // Create dataspace and dataset (use "/" prefix for root group)
      std::string fullName = "/" + name;
      hsize_t dims[1] = {1};
      auto dataspace = H5Screate_simple(1, dims, nullptr);
      auto dataset =
          H5Dcreate2(mFile.getHandle(), fullName.c_str(), memtype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      H5Dwrite(dataset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());

      H5Dclose(dataset);
      H5Sclose(dataspace);
      H5Tclose(memtype);
    }

    /**
     * @brief Save per-rank string data to a shared dataset (parallel HDF5 safe)
     * @param str The string to save
     * @param name Dataset name (shared across all ranks)
     * @param mpiRank This rank's index
     * @param nRanks Total number of MPI ranks
     *
     * Creates a dataset of size [nRanks] where each rank writes to its element.
     * This is parallel HDF5 safe - all ranks participate in collective operations.
     */
    void savePerRank(const std::string &str, const std::string &name, int mpiRank, int nRanks)
    {
      constexpr size_t LargeStringLength = 16384;

      if (str.size() >= LargeStringLength) {
        throw StringIsTooLong("String too long for HDF5 dataset '" + name + "': " + std::to_string(str.size()) +
                              " chars (max " + std::to_string(LargeStringLength) + ")");
      }

      std::vector<char> buffer(LargeStringLength, 0);
      std::strncpy(buffer.data(), str.c_str(), LargeStringLength - 1);

      // Create custom string type
      auto memtype = H5Tcopy(H5T_C_S1);
      H5Tset_size(memtype, LargeStringLength);

      // Create dataspace with nRanks elements (collective operation)
      std::string fullName = "/" + name;
      hsize_t dims[1] = {static_cast<hsize_t>(nRanks)};
      auto filespace = H5Screate_simple(1, dims, nullptr);
      auto dataset =
          H5Dcreate2(mFile.getHandle(), fullName.c_str(), memtype, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Select hyperslab for this rank's element
      hsize_t start[1] = {static_cast<hsize_t>(mpiRank)};
      hsize_t count[1] = {1};
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, count, nullptr);

      // Memory space for single element
      hsize_t memDims[1] = {1};
      auto memspace = H5Screate_simple(1, memDims, nullptr);

      // Write with independent I/O (each rank writes its own element)
      auto plist = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
      H5Pset_dxpl_mpio(plist, H5FD_MPIO_INDEPENDENT);
#endif
      H5Dwrite(dataset, memtype, memspace, filespace, plist, buffer.data());

      H5Pclose(plist);
      H5Sclose(memspace);
      H5Dclose(dataset);
      H5Sclose(filespace);
      H5Tclose(memtype);
    }

    /**
     * @brief Save per-rank double data to a shared dataset (parallel HDF5 safe)
     * @param value The value to save
     * @param name Dataset name (shared across all ranks)
     * @param mpiRank This rank's index
     * @param nRanks Total number of MPI ranks
     */
    void savePerRank(double value, const std::string &name, int mpiRank, int nRanks)
    {
      std::string fullName = "/" + name;
      hsize_t dims[1] = {static_cast<hsize_t>(nRanks)};
      auto filespace = H5Screate_simple(1, dims, nullptr);
      auto dataset = H5Dcreate2(mFile.getHandle(), fullName.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT);

      // Select hyperslab for this rank's element
      hsize_t start[1] = {static_cast<hsize_t>(mpiRank)};
      hsize_t count[1] = {1};
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, count, nullptr);

      // Memory space for single element
      hsize_t memDims[1] = {1};
      auto memspace = H5Screate_simple(1, memDims, nullptr);

      // Write with independent I/O
      auto plist = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
      H5Pset_dxpl_mpio(plist, H5FD_MPIO_INDEPENDENT);
#endif
      H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, plist, &value);

      H5Pclose(plist);
      H5Sclose(memspace);
      H5Dclose(dataset);
      H5Sclose(filespace);
    }

    /**
     * @brief Save RNG state as binary uint64_t array (compact storage)
     * @param textState String from RandomUniform::saveState() or RandomGaussian::saveState()
     * @param name Dataset name
     * @param mpiRank This rank's index
     * @param nRanks Total number of MPI ranks
     */
    void saveRNGStateBinary(const std::string &textState, const std::string &name, int mpiRank, int nRanks)
    {
      std::vector<uint64_t> binaryState;
      std::istringstream iss(textState);

      uint64_t val;
      std::string line;
      std::getline(iss, line);
      std::istringstream lineStream(line);
      while (lineStream >> val) {
        binaryState.push_back(val);
      }

      if (std::getline(iss, line) && !line.empty()) {
        std::istringstream gaussianStream(line);
        uint64_t gaussianCounter, hasCache;
        gaussianStream >> gaussianCounter >> hasCache;
        binaryState.push_back(gaussianCounter);
        binaryState.push_back(hasCache);

        if (hasCache) {
          double cachedValue;
          gaussianStream >> cachedValue;
          uint64_t cachedBits;
          std::memcpy(&cachedBits, &cachedValue, sizeof(double));
          binaryState.push_back(cachedBits);
        } else {
          binaryState.push_back(0);
        }
      }

      std::string fullName = "/" + name;
      hsize_t dims[2] = {static_cast<hsize_t>(nRanks), static_cast<hsize_t>(binaryState.size())};
      auto filespace = H5Screate_simple(2, dims, nullptr);
      auto dataset = H5Dcreate2(mFile.getHandle(), fullName.c_str(), H5T_NATIVE_UINT64, filespace, H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT);

      hsize_t start[2] = {static_cast<hsize_t>(mpiRank), 0};
      hsize_t count[2] = {1, static_cast<hsize_t>(binaryState.size())};
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, count, nullptr);

      hsize_t memDims[1] = {static_cast<hsize_t>(binaryState.size())};
      auto memspace = H5Screate_simple(1, memDims, nullptr);

      auto plist = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
      H5Pset_dxpl_mpio(plist, H5FD_MPIO_INDEPENDENT);
#endif
      H5Dwrite(dataset, H5T_NATIVE_UINT64, memspace, filespace, plist, binaryState.data());

      H5Pclose(plist);
      H5Sclose(memspace);
      H5Dclose(dataset);
      H5Sclose(filespace);
    }

    HDF5Group createOrOpenGroup(const std::string &name) { return mFile.createOrOpenGroup(name); }

    /** @brief Save the block of the dataset this rank owns, with ONE H5Dwrite per chunk.
     *
     *  The coordinates a rank contributes map to consecutive dataset indices in every dimension
     *  (see blockgeometry.h), so its whole contribution is a single stride-1 hyperslab. We stage
     *  it into one contiguous buffer and hand that to HDF5 in one call, instead of walking the
     *  first NDim-1 dimensions and writing one rod at a time: a 64^3 snapshot with 19 datasets
     *  used to cost 19*64*64 = 77824 independent 512-byte MPI-IO writes.
     *
     *  Public because nvcc's extended-lambda restrictions do not allow the DEVICE_LAMBDA below to
     *  live in a private member function.
     */
    template <typename R> void saveBlock(R r)
    {
      auto toolBox = r.getToolBox();
      constexpr size_t NDim = std::decay_t<decltype(*toolBox)>::NDim;
      using vType = typename GetGetReturnType<R>::type;

      const auto geo = computeSaveBlockGeometry<NDim>(toolBox->mLayouts.getConfigSpaceStarts(), // Local mpi offset.
                                                      toolBox->mLayouts.getConfigSpaceSizes(),  // Local mpi sizes.
                                                      toolBox->mNGridPointsVec, sparsesave, mdown, mup, mstep,
                                                      sparsesizes);
      // Nothing of this rank's subdomain is inside the saved window. The transfer mode is
      // H5FD_MPIO_INDEPENDENT, so skipping the write entirely is legal -- and is exactly what the
      // old recursion did, by never reaching writeSlices.
      if (geo.isEmpty()) return;

      const auto mLayout = toolBox->mLayouts.getConfigSpaceLayout();

      // Where the first saved coordinate of each dimension lives in local memory. The mapping is
      // coord - localStart + nGhosts, affine with slope 1, so the i-th sample along dimension d
      // sits at memFirst[d] + i * step[d].
      device::IdxArray<NDim> memFirst = geo.first;
      device::apply([&](auto... idx) { mLayout.putMemoryIndexFromSpatialLocationInto(memFirst, idx...); }, memFirst);

      // Row-major strides of the staging buffer, matching the C ordering H5Screate_simple assumes.
      // Only the dimension-0 extent is chunked below, and no stride depends on it.
      device::IdxArray<NDim> bstride{};
      bstride[NDim - 1] = 1;
      for (size_t d = NDim - 1; d > 0; --d)
        bstride[d - 1] = bstride[d] * geo.count[d];
      const device::Idx sliceElems = (NDim > 1) ? bstride[0] : 1;

      // Chunk dimension 0 so that peak staging memory stays bounded: one field's worth on device
      // plus one on host, on top of everything already resident, is an OOM waiting to happen at
      // production volumes. Every chunk is still one contiguous index box, so the file cannot tell
      // how many there were. At 64^3 doubles the loop body runs exactly once.
      const device::Idx maxSlices =
          std::max<device::Idx>(1, (device::Idx)(mStagingBudgetBytes / std::max<size_t>(1, sliceElems * sizeof(vType))));

      const device::IdxArray<NDim> starts{}; // foreach's third argument is an EXTENT, not a stop.
      const device::IdxArray<NDim> stepMem = geo.step;

      for (device::Idx c0 = 0; c0 < geo.count[0]; c0 += maxSlices) {
        const device::Idx n0 = std::min(maxSlices, geo.count[0] - c0);

        device::IdxArray<NDim> extents = geo.count;
        device::IdxArray<NDim> offsets = geo.offset;
        device::IdxArray<NDim> memBase = memFirst;
        extents[0] = n0;
        offsets[0] += c0;
        memBase[0] += c0 * geo.step[0];

        const device::Idx nElems = n0 * sliceElems;

        // Rank 1 keeps the staging view contiguous, so copyDeviceToHost takes its fast path: one
        // deep_copy, no temporary device allocation.
        device::memory::NDView<vType, 1> device_buf("SaveBlockBuffer", nElems);
        auto functor = DEVICE_LAMBDA(device::IdxArray<NDim> idx)
        {
          device::IdxArray<NDim> mem{};
          device::Idx linear = 0;
          for (size_t d = 0; d < NDim; ++d) {
            mem[d] = memBase[d] + idx[d] * stepMem[d];
            linear += idx[d] * bstride[d];
          }
          device::apply([&](const auto &...i) { device_buf(linear) = DoEval::eval(r, i...); }, mem);
        };
        device::iteration::foreach<NDim>("SaveBlockBufferFilling", starts, extents, functor);

        std::vector<vType> host(nElems);
        device::memory::copyDeviceToHost(device_buf, host.data());

        mDataset.writeSlices(host, extents, offsets);
      }
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */

    /** @brief Reject sparse limits that describe no dataset, before one is created for them.
     *
     *  setLimits cannot do this itself: it never sees the lattice, so it cannot tell whether `up`
     *  is inside it. Both rejected cases used to fail silently and late:
     *
     *  - up > nGrid made the old rod branch read `up - down` elements starting at the first owned
     *    point, i.e. straight through the ghost cells and out of the allocation, and write whatever
     *    it found to the file.
     *  - down >= up (or step < 1) gives a dataset with an empty or negative extent.
     *
     *  Note that step not dividing (up-down) is NOT rejected: it is a legitimate request, the
     *  dataset extent simply floors, and saveBlock drops the trailing coordinate that would not
     *  fit -- which is the same coordinate the old recursion failed to write anyway, only without
     *  the HDF5 error stack it printed on the way.
     */
    template <typename C> void checkSparseLimits(const C &nGrid, size_t NDim) const
    {
      if (!sparsesave) return;

      if (mdown.size() < NDim || mup.size() < NDim || mstep.size() < NDim)
        throw InvalidSparseLimits("FileSaverHDF5::setLimits was called with fewer dimensions (" +
                                  std::to_string(mup.size()) + ") than the lattice has (" + std::to_string(NDim) +
                                  ").");

      for (size_t d = 0; d < NDim; ++d) {
        const std::string where = " in dimension " + std::to_string(d) + ": down=" + std::to_string(mdown[d]) +
                                  ", up=" + std::to_string(mup[d]) + ", step=" + std::to_string(mstep[d]) +
                                  ", nGrid=" + std::to_string((device::Idx)nGrid[d]) + ".";
        if (mdown[d] < 0) throw InvalidSparseLimits("FileSaverHDF5: the lower save limit is negative" + where);
        if (mup[d] > (device::Idx)nGrid[d])
          throw InvalidSparseLimits("FileSaverHDF5: the upper save limit is beyond the lattice" + where);
        if (mdown[d] >= mup[d])
          throw InvalidSparseLimits("FileSaverHDF5: the save limits select nothing, down must be below up" + where);
        if (mstep[d] < 1) throw InvalidSparseLimits("FileSaverHDF5: the save step must be at least 1" + where);
      }
    }

    HDF5File mFile;
    HDF5Dataset mDataset;
    std::vector<device::Idx> mup, mdown, mstep;
    std::vector<device::Idx> sparsesizes;
    bool sparsesave = false;
    size_t mStagingBudgetBytes = 256ull << 20;
  };
} // namespace TempLat

#endif // HAVE_HDF5

#endif
