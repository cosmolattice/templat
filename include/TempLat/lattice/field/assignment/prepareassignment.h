#ifndef TEMPLAT_LATTICE_FIELD_ASSIGNMENT_PREPAREASSIGNMENT_H
#define TEMPLAT_LATTICE_FIELD_ASSIGNMENT_PREPAREASSIGNMENT_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include <span>
#include <utility>

#include "TempLat/lattice/algebra/helpers/confirmspace.h"
#include "TempLat/lattice/algebra/helpers/ghostreach.h"
#include "TempLat/lattice/algebra/helpers/ghostshunter.h"
#include "TempLat/lattice/memory/memorylayouts/layoutstruct.h"
#include "TempLat/lattice/memory/memorymanager.h"

namespace TempLat
{
  /** @brief Prepares the right-hand side g of an assignment that writes the fields of targets in one kernel:
   *  brings the targets and g to configuration space and walks g once to update its ghosts. Returns the range of
   *  offsets (including 0) at which g reads any target; if non-zero, assigning g in place would read already
   *  overwritten values.
   *
   *  Unit test: ctest -R test-assignmentaliasing
   **/
  template <typename T, size_t NDim, typename R>
  GhostReach prepareAssignment(R &&g, std::span<MemoryManager<T, NDim> *const> targets,
                               const LayoutStruct<NDim> &layout)
  {
    for (auto *m : targets)
      m->confirmConfigSpace(); // here the FFT may be fired
    ConfirmSpace::apply(g, layout, SpaceStateType::Configuration);

    GhostReach targetReach;
    {
      typename MemoryManager<T, NDim>::AssignmentTargetGuard guard(targets);
      GhostsHunter::apply(g);
      targetReach = guard.targetReach();
    }
    for (auto *m : targets)
      m->flagHostMirrorOutdated();
    return targetReach;
  }

  /** @brief As prepareAssignment, for a collection assigning all members in one kernel (see ComponentList), member m
   *  writing targets [m * fieldsPerMember, (m + 1) * fieldsPerMember). Walks each member separately, and also
   *  returns whether some member reads a target of another member: a kernel storing member by member would then
   *  read overwritten values.
   *
   *  Unit test: ctest -R test-sweepassignment
   **/
  template <typename T, size_t NDim, typename R>
  std::pair<GhostReach, bool> prepareCollectionAssignment(R &g, std::span<MemoryManager<T, NDim> *const> targets,
                                                          size_t fieldsPerMember, const LayoutStruct<NDim> &layout)
  {
    for (auto *m : targets)
      m->confirmConfigSpace();
    ConfirmSpace::apply(g, layout, SpaceStateType::Configuration);

    GhostReach targetReach;
    bool readsOtherMember = false;
    size_t member = 0;
    g.forEachMember([&](auto &r) {
      typename MemoryManager<T, NDim>::AssignmentTargetGuard guard(targets);
      GhostsHunter::apply(r);
      targetReach |= guard.targetReach();
      for (size_t i = 0; i < targets.size(); ++i)
        readsOtherMember = readsOtherMember || (i / fieldsPerMember != member && guard.wasRead(i));
      ++member;
    });
    for (auto *m : targets)
      m->flagHostMirrorOutdated();
    return {targetReach, readsOtherMember};
  }
} // namespace TempLat

#endif
