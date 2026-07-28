#ifndef TEMPLAT_LATTICE_IO_FILEIO_H
#define TEMPLAT_LATTICE_IO_FILEIO_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/util/exception.h"
#include "TempLat/util/rangeiteration/for_in_range.h"
#include "TempLat/util/tuplemaker.h"

#include "TempLat/lattice/IO/HDF5/filesaverhdf5.h"
#include "TempLat/lattice/IO/HDF5/fileloaderhdf5.h"

namespace TempLat
{
  MakeException(FileIOException);

  /** @brief A class which wraps HDF5 file saving and loading.
   *
   *
   *
   * Unit test: ctest -R test-fileio
   **/
  template <size_t NDim = 3> class FileIO
  {
  public:
    // Put public methods here. These should change very little over time.
    FileIO() = default;

#if defined(HAVE_HDF5) && defined(HAVE_MPI)
    /** @brief Point HDF5 collective I/O at the communicator the lattice is decomposed over.
     *
     *  Defaults to MPI_COMM_WORLD. Pass `toolBox->getMPIGroup().getBaseComm()` (or whatever
     *  communicator the lattice was decomposed over) when running on a sub-communicator —
     *  otherwise the collective barriers inside HDF5File deadlock against ranks that never
     *  participate, and per-rank dataset offsets index the wrong ranks. Call before open/create.
     */
    void setComm(MPI_Comm comm)
    {
      saver.setComm(comm);
      loader.setComm(comm);
    }
#endif

    template <class R>
      requires HasStaticGetter<typename std::decay_t<R>>
    void save(R &&r)
    {
#ifdef HAVE_HDF5
      using nakedR = std::decay_t<R>;
      for_in_range<number_to_skip_as_tuple<nakedR>::value, nakedR::size>(
          [&](auto i) { save(nakedR::Getter::get(r, i)); });
#else
      throw(FileIOException("You tried to save an object to a file, but the HDF5 library is not available. Make sure "
                            "you have it installed and that you compiled CosmoLattice with it."));
#endif
    }

    template <class R>
      requires(!HasStaticGetter<typename std::decay_t<R>>)
    void save(R &&r)
    {
#ifdef HAVE_HDF5
      saver.save(r);
#else
      throw(FileIOException("You tried to save an object to a file, but the HDF5 library is not available. Make sure "
                            "you have it installed and that you compiled CosmoLattice with it."));
#endif
    }

    template <class R>
      requires HasStaticGetter<typename std::decay_t<R>>
    void load(R &&r)
    {
#ifdef HAVE_HDF5
      using nakedR = std::decay_t<R>;
      for_in_range<number_to_skip_as_tuple<nakedR>::value, nakedR::size>(
          [&](auto i) { load(nakedR::Getter::get(r, i)); });
#else
      throw(FileIOException("You tried to load an object from a file, but the HDF5 library is not available. Make sure "
                            "you have it installed and that you compiled CosmoLattice with it."));
#endif
    }

    template <class R>
      requires(!HasStaticGetter<typename std::decay_t<R>>)
    void load(R &&r)
    {
#ifdef HAVE_HDF5
      loader.load(r);
#else
      throw(FileIOException("You tried to load an object from a file, but the HDF5 library is not available. Make sure "
                            "you have it installed and that you compiled CosmoLattice with it."));
#endif
    }

    template <typename R> auto setSaverLimits(R down, R up, R step)
    {
#ifdef HAVE_HDF5
      saver.setLimits(down, up, step, NDim);
#else
      throw(FileIOException(
          "You tried to set a subvolume for the HDF5 snapshots, but the HDF5 library is not available. Make sure "
          "you have it installed and that you compiled CosmoLattice with it."));
#endif
    }

#ifdef HAVE_HDF5
    FileSaverHDF5 saver;
    FileLoaderHDF5 loader;
#endif
  };
} // namespace TempLat

#endif
