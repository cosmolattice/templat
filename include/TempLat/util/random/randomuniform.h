#ifndef TEMPLAT_UTIL_RANDOM_RANDOMUNIFORM_H
#define TEMPLAT_UTIL_RANDOM_RANDOMUNIFORM_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Franz R. Sattler,  Year: 2025

#include <cstdint>
#include <random>
#include <sstream>

#include "TempLat/util/hash/keccakhash.h"
#include "TempLat/parallel/device.h"
#include "TempLat/parallel/device_memory.h"

#include <Random123/philox.h>

namespace TempLat
{
  /** @brief A class which gives pseudo random counter-based rng, based on a string random seed, stable across
   * platforms.
   *
   *
   * Unit test: ctest -R test-randomuniform
   **/
  template <typename T, typename RNG = r123::Philox2x64> class RandomUniform
  {
    using INT = typename RNG::ctr_type::value_type;

  public:
    // Put public methods here. These should change very little over time.

    using IntegerType = INT;

    RandomUniform(const std::string &stringSeed)
        : mStringSeed(stringSeed), mHashSeed(KeccakHash::compute(stringSeed)),
          mSeed(static_cast<INT>((uint64_t)mHashSeed))
    {
    }

    /**
     * @brief Serializes the complete RNG state to a string
     * @return String containing serialized mt19937_64 state and counter
     */
    std::string saveState() const
    {
      std::ostringstream oss;
      oss << *mStringSeed;
      return oss.str();
    }

    /**
     * @brief Restores RNG state from a serialized string
     * @param state String produced by saveState()
     */
    void loadState(const std::string &state)
    {
      std::istringstream iss(state);
      iss >> *mStringSeed;
      mHashSeed = KeccakHash::compute(*mStringSeed);
      mSeed = static_cast<INT>((uint64_t)mHashSeed);
    }

    auto getSeed() const { return mSeed; }

    const std::string &getSeedString() const { return *mStringSeed; }

    DEVICE_FORCEINLINE_FUNCTION
    T get(INT r, INT c, INT gen) const { return getPair(r, c, gen)[0]; }

    DEVICE_FORCEINLINE_FUNCTION
    device::array<T, 2> getPair(INT r, INT c, INT gen) const
    {
      const RNG rng;

      // create a counter and key for the generator
      const typename RNG::ctr_type counter = {{r, c}};
      const typename RNG::key_type key = {{mSeed + gen}};
      // draw a pair of numbers
      const typename RNG::ctr_type result = rng(counter, key);

      return {{integer_to_float(result[0]), integer_to_float(result[1])}};
    }

    template <typename R> DEVICE_FORCEINLINE_FUNCTION T integer_to_float(R value) const
    {
      return (static_cast<T>(value) - static_cast<T>(std::numeric_limits<R>::min())) /
             static_cast<T>(std::numeric_limits<R>::max());
    }

    /** @brief For testing purposes, we need to compare prng's, specifically their seeds. */
    friend bool operator==(const RandomUniform &a, const RandomUniform &b) { return a.getSeed() == b.getSeed(); }

    friend std::ostream &operator<<(std::ostream &ostream, const RandomUniform &pr)
    {
      ostream << "RandomUniform - seed string: \"" << pr.getSeedString() << "\" - seed value: " << pr.getSeed() << "\n";
      return ostream;
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    device::memory::host_string mStringSeed;
    KeccakHash::ResultType mHashSeed;
    INT mSeed;
  };
} // namespace TempLat

#endif
