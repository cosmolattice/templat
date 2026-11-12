/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

/* Pure-gauge evolution of SU2_NGROUPS independent SU(2) gauge fields, stored in a VectorFieldCollection:
 * does assigning all groups in one kernel pay off?
 *
 *   default            one kernel per group (VectorField assignment)
 *   -DGROUPS_FUSED     one kernel for all groups and directions, emulating a fused VectorFieldCollection assignment
 *   -DGROUPS_PLAQGRAD  one kernel per group with plaqGrad, for reference
 *   -DSU2_NGRID, -DSU2_NSTEPS, -DSU2_NGROUPS set the lattice size, step count and number of groups.
 *
 * The kick is the one of su2_evolution.cpp, the drift U_i = exp(dt Pi_i) U_i.
 */

#include "TempLat/session/sessionguard.h"
#include "TempLat/util/benchmark.h"

#include "TempLat/lattice/algebra/su2algebra/su2algebra.h"
#include "TempLat/lattice/algebra/gaugealgebra/plaquette.h"
#include "TempLat/lattice/algebra/gaugealgebra/plaquetteback.h"
#include "TempLat/lattice/algebra/gaugealgebra/plaquettegradient.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/listoperators/listoperators.h"
#include "TempLat/lattice/field/collections/vectorfieldcollection.h"
#include "TempLat/lattice/field/assignment/assignviews.h"
#include "TempLat/lattice/field/assignment/componentlist.h"
#include "TempLat/lattice/field/assignment/prepareassignment.h"
#include "TempLat/lattice/measuringtools/averager.h"
#include "TempLat/util/rangeiteration/for_in_range.h"
#include "TempLat/util/rangeiteration/make_list_tag.h"
#include "TempLat/util/rangeiteration/sum_in_range.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <utility>

#ifndef SU2_NGRID
#define SU2_NGRID 64
#endif
#ifndef SU2_NSTEPS
#define SU2_NSTEPS 100
#endif
#ifndef SU2_NGROUPS
#define SU2_NGROUPS 2
#endif

namespace
{
  using namespace TempLat;

  template <typename A> struct FieldCount;
  template <typename F, size_t K> struct FieldCount<device::array<F, K>> {
    static constexpr size_t value = K;
  };

  /** @brief Assigns rs[m] to targets[m], all in one kernel, as FieldCollection::operator= does for its members. */
  template <typename T, size_t NDim, typename F, size_t M, typename... Rs>
  void assignFused(const std::array<F, M> &targets, const Rs &...rs)
  {
    static_assert(sizeof...(Rs) == M);
    ComponentList<T, Rs...> g(rs...);
    constexpr size_t K = FieldCount<std::remove_cvref_t<decltype(targets[0].assignedFields())>>::value;
    const auto first = targets[0].assignedFields()[0];
    device::array<decltype(first.getView()), M * K> views;
    std::array<MemoryManager<T, NDim> *, M * K> managers;
    for (size_t m = 0; m < M; ++m)
      for (size_t k = 0; k < K; ++k) {
        const auto f = targets[m].assignedFields()[k];
        views[m * K + k] = f.getView();
        managers[m * K + k] = f.getMemoryManager().get();
      }
    const auto &layout = first.getLayout();
    const auto [reach, readsOtherMember] = prepareCollectionAssignment<T, NDim>(g, managers, K, layout);
    assignViews<F::firstAssignedComponent>("GroupsAssign", views, g, layout, reach, first.toolBox()->mSweepBuffer,
                                           readsOtherMember);
    for (auto *m : managers)
      m->setGhostsAreStale();
  }
} // namespace

int main(int argc, char **argv)
{
  using namespace TempLat;

  SessionGuard guard(argc, argv, false);

  constexpr size_t NDim = 3;
  using T = double;
  constexpr size_t nGrid = SU2_NGRID;
  constexpr size_t nSteps = SU2_NSTEPS;
  constexpr int G = SU2_NGROUPS;
  constexpr T dt = 0.01;
  constexpr T amp = 0.5;

  auto toolBox = MemoryToolBox<NDim>::makeShared(nGrid, 1, false);
  toolBox->unsetVerbose();

  VectorFieldCollection<SU2Field<T, NDim>, G> U("U", toolBox);
  VectorFieldCollection<SU2LieAlgebraField<T, NDim>, G> Pi("Pi", toolBox);

  SpatialCoordinate<NDim> x(toolBox);
  const T w = 2.0 * M_PI / nGrid;
  for_in_range<0, G>([&](auto g) {
    const T a = amp * (1 + 0.2 * decltype(g)::value);
    U(g)(1_c) = MakeSU2(c, IfElse(c == 0_c, cos(a * sin(w * x(1_c))), IfElse(c == 1_c, sin(a * sin(w * x(1_c))), 0.)));
    U(g)(2_c) = MakeSU2(c, IfElse(c == 0_c, cos(a * sin(w * x(2_c))), IfElse(c == 2_c, sin(a * sin(w * x(2_c))), 0.)));
    U(g)(3_c) = MakeSU2(c, IfElse(c == 0_c, cos(a * sin(w * x(3_c))), IfElse(c == 3_c, sin(a * sin(w * x(3_c))), 0.)));
    for_in_range<1, NDim + 1>([&](auto i) {
      Pi(g)(i)(1_c) = T(0);
      Pi(g)(i)(2_c) = T(0);
      Pi(g)(i)(3_c) = T(0);
    });
  });

  auto energy = [&]() {
    T e = 0;
    for_in_range<0, G>([&](auto g) {
      e += 2.0 * Total(i, 1, NDim, average(pow<2>(Pi(g)(i)(1_c)) + pow<2>(Pi(g)(i)(2_c)) + pow<2>(Pi(g)(i)(3_c))));
      e += 2.0 * Total(i, 1, NDim, Total(j, 1, NDim, IfElse(i != j, average(1.0 - plaq(U(g), i, j).SU2Get(0_c)), 0.0)));
    });
    return e;
  };
  const T e0 = energy();

  auto kick = [&](auto g, auto i) {
    return Pi(g)(i) - dt * Total(j, 1, NDim, IfElse(i != j, plaq(U(g), i, j) - plaqBack(U(g), i, j), ZeroType()));
  };
  auto drift = [&](auto g, auto i) { return exp(dt * Pi(g)(i)) * U(g)(i); };
#if defined(GROUPS_FUSED)
  // Member m of the fused assignment is group m / NDim, direction m % NDim + 1.
  auto group = [](auto m) { return Tag<int(decltype(m)::value / NDim)>(); };
  auto dir = [](auto m) { return Tag<int(decltype(m)::value % NDim + 1)>(); };
#endif

  Benchmark bench([&](Benchmark::Measurer &measurer) {
    for (size_t step = 0; step < nSteps; ++step) {
      measurer.measure("Pi kick", [&]() {
#if defined(GROUPS_FUSED)
        [&]<size_t... m>(std::index_sequence<m...>) {
          assignFused<T, NDim>(std::array{Pi(group(Tag<m>()))(dir(Tag<m>()))...},
                               kick(group(Tag<m>()), dir(Tag<m>()))...);
        }(std::make_index_sequence<G * NDim>());
#elif defined(GROUPS_PLAQGRAD)
        for_in_range<0, G>([&](auto g) { Pi(g) = Pi(g) - dt * plaqGrad(U(g)); });
#else
        for_in_range<0, G>([&](auto g) { Pi(g) = MakeVector(i, 1, NDim, kick(g, i)); });
#endif
        device::iteration::fence();
      });
      measurer.measure("U drift", [&]() {
#if defined(GROUPS_FUSED)
        [&]<size_t... m>(std::index_sequence<m...>) {
          assignFused<T, NDim>(std::array{U(group(Tag<m>()))(dir(Tag<m>()))...},
                               drift(group(Tag<m>()), dir(Tag<m>()))...);
        }(std::make_index_sequence<G * NDim>());
#else
        for_in_range<0, G>([&](auto g) { U(g) = MakeVector(i, 1, NDim, drift(g, i)); });
#endif
        device::iteration::fence();
      });
    }
  });
  bench.run(1);
  bench.print();
  bench.log("su2_groups");

  const T e1 = energy();
  std::printf("energy: initial = %.12f  final = %.12f  relative drift = %.3e\n", e0, e1, std::abs(e1 / e0 - 1.0));
}
