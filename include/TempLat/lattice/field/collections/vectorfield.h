#ifndef COSMOINTERFACE_3D_VECTORFIELD3D_H
#define COSMOINTERFACE_3D_VECTORFIELD3D_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Adrien Florio,  Year: 2019

#include "TempLat/lattice/field/collections/fieldcollection.h"
#include "TempLat/lattice/algebra/helpers/getndim.h"

namespace TempLat
{
  /** @brief A class which
   * Field collections. Allows to have vector fields, index starting from one.
   *
   *
   * Unit test: ctest -R test-vectorfield3d
   *
   * @vocab-summary A vector-valued field — one component per spatial dimension, indexed from 1. The natural
   * type for a gauge link $U_\mu$.
   * @vocab-signature VectorField<FieldType> Us("Us", toolBox);
   * @vocab-tags Collection
   **/
  template <class Arg, bool flatAssign = false>
  using VectorField = FieldCollection<Arg, GetNDim::get<Arg>(), flatAssign, 1>;
} // namespace TempLat

#endif
