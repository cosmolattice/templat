
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2019
#include "TempLat/lattice/algebra/helpers/getvectorcomponent.h"
#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/algebra/spacestateinterface.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/lattice/field/collections/fieldcollection.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/helpers/haseval.h"
#include "TempLat/lattice/algebra/coordinates/wavenumber.h"

namespace TempLat
{
  template <int N, typename T> struct GetVectorComponentHelperTester {
    static void Test(TDDAssertion &tdd);
  };

  template <int N, typename T> inline void GetVectorComponentHelperTester<N, T>::Test(TDDAssertion &tdd)
  {
    auto toolBox = MemoryToolBox<3>::makeShared(32, 1);
    toolBox->setVerbose();
    FieldCollection<Field<T, 3>, 3> fc("abcdefg", toolBox);

    fc[1_c] = 1;
    fc[2_c].inFourierSpace() = 2;
    tdd.verify(fc[1_c].isFourierSpace() == false);
    tdd.verify(fc[2_c].isFourierSpace() == true);

    auto test1 = GetVectorComponentHelper<1, FieldCollection<Field<T, 3>, 3>>(fc);
    auto test2 = GetVectorComponentHelper<2, FieldCollection<Field<T, 3>, 3>>(fc);

    auto test_wavenumber = GetVectorComponentHelper<1, WaveNumber<3>>(WaveNumber<3>(toolBox));

    tdd.verify(HasEvalMethod<decltype(test_wavenumber)>);
    static_assert(HasEvalMethod<decltype(test_wavenumber)>);

    test1.confirmSpace(toolBox->mLayouts.getConfigSpaceLayout(), SpaceStateType::Configuration);
    tdd.verify(fc[2_c].isFourierSpace() == true);
    test2.confirmSpace(toolBox->mLayouts.getConfigSpaceLayout(), SpaceStateType::Configuration);
    tdd.verify(fc[2_c].isFourierSpace() == false);

    fc[1_c].setGhostsAreStale();
    test2.confirmGhostsUpToDate();
    tdd.verify(fc[1_c].areGhostsStale() == true);
    test1.confirmGhostsUpToDate();
    tdd.verify(fc[1_c].areGhostsStale() == false);

    tdd.verify(test1.toString() == "abcdefg_1(x)");
  }
} // namespace TempLat

namespace
{
  TempLat::TDDContainer<TempLat::GetVectorComponentHelperTester<1, double>> test;
#ifdef HAVE_FFTFLOAT
  TempLat::TDDContainer<TempLat::GetVectorComponentHelperTester<1, float>> test2;
#endif
}
