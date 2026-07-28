
/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Wessel Valkenburg, Year: 2019
#include "TempLat/lattice/memory/memorylayouts/hermitianpartners.h"
#include "TempLat/util/tdd/tdd.h"

namespace TempLat
{

  template <size_t NDim> struct HermitianPartnersTester {
    static void Test(TDDAssertion &tdd);
  };

  template <size_t NDim> void HermitianPartnersTester<NDim>::Test(TDDAssertion &tdd)
  {
    HermitianPartners<4> hp({16, 16, 16, 16});
    /* simply test the default implementation */
    tdd.verify(hp.qualify(device::IdxArray<4>{{0, 0, 0, 0}}) == HermitianRedundancy::none);

    device::IdxArray<4> global{{16, 16, 16, 16}}, hpvec{{0, 0, 0, 0}};
    hp.putHermitianPartner(global, hpvec);
    bool correct = true;
    for (device::Idx i = 0; i < 4; ++i)
      correct = correct && global[i] == hpvec[i];
    tdd.verify(correct);
  }

} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::HermitianPartnersTester<4>> test;
}
