#ifndef TEMPLAT_LATTICE_FIELD_ASSIGNMENT_PREPAREASSIGNMENT_H
#define TEMPLAT_LATTICE_FIELD_ASSIGNMENT_PREPAREASSIGNMENT_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Franz R. Sattler, Year: 2026

#include <span>

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
} // namespace TempLat

#endif
