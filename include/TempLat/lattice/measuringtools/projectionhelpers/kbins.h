#ifndef TEMPLAT_LATTICE_MEASUREMENTS_PROJECTIONHELPERS_HALKKBINS_H
#define TEMPLAT_LATTICE_MEASUREMENTS_PROJECTIONHELPERS_HALKKBINS_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio, Year: 2019

#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionresult.h"

namespace TempLat
{
  /**
   *  Simple convenience class which maps any bining to one integer ks. Useful for test purposes.
   *
   *
   * Unit test: ctest -R test-kbins
   **/
  template <typename T> class KBins
  {
  public:
    // Put public methods here. These should change very little over time.
    KBins(const device::Idx &nGrid, const RadialProjectionResult<T> &res, const T &pKir)
        : ks((int)(std::pow(3, 0.5) / 2.0 * (T)nGrid)), vs((int)(std::pow(3, 0.5) / 2.0 * (T)nGrid)),
          ms((int)(std::pow(3, 0.5) / 2.0 * (T)nGrid)), kir(pKir)
    {
      size_t ind = 0;
      for (auto &&it : res) {
        if (!AlmostEqual(it.getBin().average, 0)) {
          // Clamp both ends: round(avg)-1 can be negative for a first-bin mean in (0, 0.5), which would
          // underflow the size_t index and write out of bounds.
          ind = (size_t)std::max(0.0, std::min(std::round((it.getBin().average)) - 1, ks.size() - 1.0));
          //  say<<ind<<"\n";
          vs[ind] += it.getValue().average * it.getValue().multiplicity;
          ms[ind] += it.getValue().multiplicity;
        }
      }
      for (size_t i = 0; i < ks.size(); ++i) {
        ks[i] = i + 1;
        if (ms[i] != 0) vs[i] /= ms[i]; // leave empty bins at zero instead of dividing by zero
      }
    }

    friend std::ostream &operator<<(std::ostream &ostream, const KBins &rpr)
    {
      // ostream << "#ks values multiplicities" << "\n";
      for (size_t i = 0; i < rpr.ks.size(); ++i) {
        ostream << rpr.ks[i] * rpr.kir << " " << rpr.vs[i] << " " << rpr.ms[i] << "\n";
      }
      // ostream <<"\n";
      return ostream;
    }

  private:
    /* Put all member variables and private methods here. These may change arbitrarily. */
    std::vector<T> ks; // ks
    std::vector<T> vs; // spectrum or whatever is binned
    std::vector<T> ms; // multiplicities, to check
    T kir;
  };
} // namespace TempLat

#endif
