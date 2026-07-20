#ifndef TEMPLAT_PARALLEL_MPI_COMM_MPICARTESIANGROUP_H
#define TEMPLAT_PARALLEL_MPI_COMM_MPICARTESIANGROUP_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg,  Year: 2019

#include "TempLat/util/exception.h"

#include "TempLat/parallel/mpi/comm/mpicommreference.h"
#include "TempLat/parallel/mpi/comm/mpidomainsplit.h"

namespace TempLat
{

  MakeException(MPICartesianGroupException);

  /** @brief A class keeps the books regarding the process layout
   * relative to the physics lattice.
   *
   * Unit test: ctest -R test-mpicartesiangroup
   **/
  class MPICartesianGroup
  {
  public:
    // Put public methods here. These should change very little over time.
    /** @brief Constructor that adopts an externally-determined topology, normally the FFT
     *        backend's (see `FFTTopology`, built by `FFTMPIDomainSplit::makeMPIGroup`).
     *
     *  The Cartesian communicator is created over `topologyGroup` with **no reordering**, so this
     *  group's rank->coordinate mapping is exactly that communicator's. `expectedCoords`, when
     *  non-empty, is then checked against the coordinates MPI actually hands back, and a mismatch
     *  throws.
     *
     *  This matters because the local starts installed into the layout come from the backend's
     *  coordinates, while ghost neighbours come from `MPI_Cart_shift` on the communicator built
     *  here. These used to be two independently-reordered topologies reconciled only by comparing
     *  grid *shape*; when they disagreed, subdomain boundaries exchanged the wrong data silently.
     *
     *  \param topologyGroup the MPI Comm whose rank order defines the topology.
     *  \param nDimensions the dimensionality of the Cartesian setup.
     *  \param decomposition the grid shape (rank counts per lattice dimension).
     *  \param expectedCoords this rank's coordinates as the topology's owner reports them, or
     *         empty to skip the cross-check.
     *  \param (optional)periodic per-dimension periodicity (1=true, 0=false). Default all periodic.
     */
    MPICartesianGroup(MPICommReference topologyGroup, ptrdiff_t nDimensions, std::vector<int> decomposition,
                      std::vector<int> expectedCoords = std::vector<int>(),
                      std::vector<int> periodic = std::vector<int>())
        : mBaseGroup(topologyGroup), mCartesianGroup(mBaseGroup), mNDimensions(nDimensions),
          mDecomposition(decomposition), mPeriodic(periodic), mExpectedCoords(expectedCoords)
    {
      verifyInput();

      createGroups();

      verifyCoordsMatchTopology();
    }

    /** @brief convenience short-hand constructor with MPI_COMM_WORLD */
    MPICartesianGroup(ptrdiff_t nDimensions, std::vector<int> decomposition)
        : MPICartesianGroup(MPICommReference(), nDimensions, decomposition)
    {
    }

    /** @brief Get the MPI_Comm already. */
    MPI_Comm getComm() const { return mCartesianGroup; }
    MPICommReference getBaseComm() const { return mBaseGroup; }

    int getRank() { return mCartesianGroup.getRank(); }

    ptrdiff_t getNDimensions() const { return mNDimensions; }

    /** @brief Returns the position of the current process in the cartesian grid. Values are relative
     *         to the rank layout, nothing else. */
    const std::vector<int> &getPosition() const { return mSelfPosition; }

    /** @brief The decomposition (grid shape) of this group. */
    const std::vector<int> &getDecomposition() const { return mDecomposition; }

    auto size() const { return mBaseGroup.size(); }

    friend std::ostream &operator<<(std::ostream &ostream, MPICartesianGroup mGr)
    {
      ostream << "MPICartesianGroup:\n  decomposition: " << mGr.getDecomposition() << "\n  this rank: " << mGr.getRank()
              << "\n  this position: " << mGr.getPosition() << "\n";
      return ostream;
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */

    /**
     * @brief The BaseGroup is the MPI_Comm that is given to MPICartesianGroup in the constructor, and from which we
     * create the Cartesian groups.
     *
     */
    MPICommReference mBaseGroup;
    /**
     * @brief The CartesianGroup is the MPI_Comm that is created from the BaseGroup, and is split up in the number of
     * dimensions that were passed in the constructor.
     */
    MPICommReference mCartesianGroup;

    const ptrdiff_t mNDimensions;
    std::vector<int> mDecomposition;
    std::vector<int> mPeriodic;
    std::vector<int> mExpectedCoords;
    std::vector<int> mSelfPosition;

    std::vector<int> fetchPosition(int ofRank)
    {
      std::vector<int> result(mDecomposition.size(), 0);
#ifdef HAVE_MPI
      // sayMPI << this << "Getting coords: " << ofRank << ", " << result << " nd: " <<mNDimensions<< "\n";
      if (0 != MPI_Cart_coords(mCartesianGroup, ofRank, (int)result.size(), result.data())) {
        throw MPICartesianGroupException("Could not get result from MPI_Cart_coords.");
      }
#endif
      return result;
    }

    void createGroups()
    {
#ifdef HAVE_MPI
      /* No reordering. The rank order of mBaseGroup is the topology, so the Cartesian
       * communicator must preserve it: the layout's local starts are derived from that order,
       * and ghost neighbours are derived from MPI_Cart_shift on the communicator built here.
       * Letting MPI permute ranks lets those two disagree. */
      mCartesianGroup = createOneGroup(mBaseGroup, mNDimensions, false);
#endif

      mSelfPosition = fetchPosition(mCartesianGroup.rank());
    }

    MPICommReference createOneGroup(MPICommReference fromGroup, ptrdiff_t nDim, bool reorder)
    {
      MPI_Comm newComm;

#ifdef HAVE_MPI
      /* Build from fromGroup, not mBaseGroup: this used to ignore its own argument. */
      MPI_Cart_create(fromGroup, nDim, mDecomposition.data(), mPeriodic.data(), reorder, &newComm);
#else
      (void)fromGroup;
      (void)nDim;
      (void)reorder;
      newComm = MPI_COMM_WORLD;
#endif
      return MPICommReference(newComm);
    }

    /** @brief Confirm the coordinates MPI computed are the ones the topology's owner expects.
     *
     *  With `reorder = false` this holds by construction, so a failure here means either that an
     *  MPI implementation permuted ranks anyway or that the supplied grid shape does not describe
     *  the supplied communicator. Both silently corrupt subdomain boundaries, so fail loudly.
     */
    void verifyCoordsMatchTopology() const
    {
      if (mExpectedCoords.empty()) return;

      if ((ptrdiff_t)mExpectedCoords.size() != mNDimensions)
        throw MPICartesianGroupException("Expected coordinates have ", mExpectedCoords.size(), " entries, but this is a ",
                                         mNDimensions, "-dimensional group.");

      for (ptrdiff_t i = 0; i < mNDimensions; ++i) {
        if (mSelfPosition[i] != mExpectedCoords[i])
          throw MPICartesianGroupException(
              "Cartesian coordinates disagree with the topology that produced this group: MPI places this rank at ",
              mSelfPosition, " but the topology's owner reports ", mExpectedCoords, " (decomposition ", mDecomposition,
              "). The layout's local starts follow the latter while ghost exchange follows the former, so continuing "
              "would silently exchange the wrong data at subdomain boundaries.");
      }
    }

    void verifyInput()
    {
      while ((ptrdiff_t)mDecomposition.size() < mNDimensions)
        mDecomposition.push_back(1);
      while (mPeriodic.size() < mDecomposition.size())
        mPeriodic.push_back(1);

      if ((ptrdiff_t)mDecomposition.size() > mNDimensions) mDecomposition.resize(mNDimensions);
      if ((ptrdiff_t)mPeriodic.size() > mNDimensions) mPeriodic.resize(mNDimensions);

      /* The decomposition must tile the communicator exactly. This used to "repair" a mismatch
       * by overwriting mDecomposition[0] and then zeroing trailing entries from the right, which
       * silently mangled an explicitly-pinned backend grid into something the backend would not
       * use — precisely the disagreement this class now exists to prevent. */
      ptrdiff_t groupSize = mBaseGroup.size();
      ptrdiff_t product = 1;
      for (auto &&it : mDecomposition)
        product *= it;

      if (product != groupSize)
        throw MPICartesianGroupException("Decomposition ", mDecomposition, " describes ", product,
                                         " ranks but the communicator has ", groupSize,
                                         ". The decomposition must tile the communicator exactly; build the group via "
                                         "FFTMPIDomainSplit::makeMPIGroup so it matches the FFT backend.");
    }
  };
} // namespace TempLat

#endif
