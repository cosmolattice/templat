#ifndef TEMPLAT_H
#define TEMPLAT_H

/* This file is part of TempLat, available at https://cosmolattice.github.io/templat .
   Copyright 2021-2026 The TempLat authors, see AUTHORS.md.
   Released under the MIT license, see LICENSE.md. */

/** \file The umbrella header of TempLat.
 *
 * TempLat is header-only, so a program can pull in the whole library with a single
 *
 *     #include <TempLat.h>
 *
 * and never has to know which of the ~400 individual headers a given class lives in.
 * Everything below is guarded so that it compiles in every supported configuration:
 * the MPI, HDF5, ParaFaFT and KokkosFFT parts disable themselves when TempLat was not
 * built against those libraries.
 *
 * Including everything is not free at compile time. Translation units that only touch a
 * small corner of the library are still free to include the individual headers directly,
 * exactly as before -- this header adds a convenience, it does not replace anything.
 *
 * Two groups of headers are deliberately *not* included here:
 *
 *  - The FFT backends under TempLat/fft/external/ and the device backends under
 *    TempLat/parallel/devices/. These are implementation details that are mutually
 *    exclusive and picked for you by TempLat/fft/fftlibraryselector.h and
 *    TempLat/parallel/device.h respectively, both of which *are* included below.
 *  - The unit-test harness under TempLat/util/tdd/. Include TempLat/util/tdd/tdd.h
 *    yourself if you are writing tests against TempLat.
 *
 * Unit test: ctest -R test-templat
 */

// ---------------------------------------------------------------------------
// Session management
// ---------------------------------------------------------------------------
#include "TempLat/session/sessionguard.h"

// ---------------------------------------------------------------------------
// Parallel backends: device abstraction, threading and MPI
// ---------------------------------------------------------------------------
#include "TempLat/parallel/device.h"
#include "TempLat/parallel/device_guard.h"
#include "TempLat/parallel/device_iteration.h"
#include "TempLat/parallel/device_memory.h"
#include "TempLat/parallel/mpi/cartesian/mpicartesianexchange.h"
#include "TempLat/parallel/mpi/cartesian/mpicartesiangroup.h"
#include "TempLat/parallel/mpi/cartesian/mpicartesianneighbours.h"
#include "TempLat/parallel/mpi/cartesian/mpicartesianneighbourssingledimension.h"
#include "TempLat/parallel/mpi/comm/exchange/mpiallreduce.h"
#include "TempLat/parallel/mpi/comm/exchange/mpisendreceive.h"
#include "TempLat/parallel/mpi/comm/mpicommreference.h"
#include "TempLat/parallel/mpi/comm/mpidomainsplit.h"
#include "TempLat/parallel/mpi/mpitags.h"
#include "TempLat/parallel/mpi/mpitypeconstants.h"
#include "TempLat/parallel/mpi/session/mpiguard.h"
#include "TempLat/parallel/threadsettings.h"

// ---------------------------------------------------------------------------
// Runtime parameters: input files and the command line
// ---------------------------------------------------------------------------
#include "TempLat/parameters/filereader.h"
#include "TempLat/parameters/multipleparametergetter.h"
#include "TempLat/parameters/pairmaker.h"
#include "TempLat/parameters/parametergetter.h"
#include "TempLat/parameters/parameterparser.h"
#include "TempLat/parameters/stringconverter.h"

// ---------------------------------------------------------------------------
// Fast Fourier transforms
// ---------------------------------------------------------------------------
#include "TempLat/fft/fftdecomposition.h"
#include "TempLat/fft/fftlibraryinterface.h"
#include "TempLat/fft/fftlibraryselector.h"
#include "TempLat/fft/fftmpidomainsplit.h"
#include "TempLat/fft/fftnormalization.h"
#include "TempLat/fft/ffttopology.h"

// ---------------------------------------------------------------------------
// Lattice parameters
// ---------------------------------------------------------------------------
#include "TempLat/lattice/latticeparameters.h"

// ---------------------------------------------------------------------------
// Lattice memory: layouts, blocks and the memory toolbox
// ---------------------------------------------------------------------------
#include "TempLat/lattice/memory/memoryblock.h"
#include "TempLat/lattice/memory/memorylayouts/checkerboardlayout.h"
#include "TempLat/lattice/memory/memorylayouts/fftlayoutstruct.h"
#include "TempLat/lattice/memory/memorylayouts/hermitianpartners.h"
#include "TempLat/lattice/memory/memorylayouts/hermitianredundancy.h"
#include "TempLat/lattice/memory/memorylayouts/hermitianvalueaccounting.h"
#include "TempLat/lattice/memory/memorylayouts/layoutstruct.h"
#include "TempLat/lattice/memory/memorylayouts/layoutstructglobal.h"
#include "TempLat/lattice/memory/memorylayouts/layoutstructlocal.h"
#include "TempLat/lattice/memory/memorylayouts/layoutstructlocaltransposed.h"
#include "TempLat/lattice/memory/memorylayouts/transpositionmap.h"
#include "TempLat/lattice/memory/memorylayoutstate.h"
#include "TempLat/lattice/memory/memorymanager.h"
#include "TempLat/lattice/memory/memorytoolbox.h"
#include "TempLat/lattice/memory/triplestatelayouts.h"
#include "TempLat/lattice/memory/verbositylevels.h"

// ---------------------------------------------------------------------------
// Ghost cells
// ---------------------------------------------------------------------------
#include "TempLat/lattice/ghostcells/ghostbuster.h"
#include "TempLat/lattice/ghostcells/ghoststatekeeper.h"
#include "TempLat/lattice/ghostcells/ghostsubarray.h"
#include "TempLat/lattice/ghostcells/ghostsubarraymap.h"
#include "TempLat/lattice/ghostcells/ghostupdater.h"

// ---------------------------------------------------------------------------
// Fields and field collections
// ---------------------------------------------------------------------------
#include "TempLat/lattice/field/abstractfield.h"
#include "TempLat/lattice/field/assignablefieldcollection.h"
#include "TempLat/lattice/field/collections/fieldcollection.h"
#include "TempLat/lattice/field/collections/vectorfield.h"
#include "TempLat/lattice/field/collections/vectorfieldcollection.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/field/helpers/flattener.h"
#include "TempLat/lattice/field/helpers/hasastuplecat.h"
#include "TempLat/lattice/field/views/fieldviewconfig.h"
#include "TempLat/lattice/field/views/fieldviewfourier.h"

// ---------------------------------------------------------------------------
// Field algebra: elementary operators
// ---------------------------------------------------------------------------
#include "TempLat/lattice/algebra/operators/absolutevalue.h"
#include "TempLat/lattice/algebra/operators/acos.h"
#include "TempLat/lattice/algebra/operators/add.h"
#include "TempLat/lattice/algebra/operators/arg.h"
#include "TempLat/lattice/algebra/operators/arg2.h"
#include "TempLat/lattice/algebra/operators/asinh.h"
#include "TempLat/lattice/algebra/operators/binaryoperator.h"
#include "TempLat/lattice/algebra/operators/complexconjugate.h"
#include "TempLat/lattice/algebra/operators/cosh.h"
#include "TempLat/lattice/algebra/operators/cosine.h"
#include "TempLat/lattice/algebra/operators/diracdeltafunction.h"
#include "TempLat/lattice/algebra/operators/divide.h"
#include "TempLat/lattice/algebra/operators/exponential.h"
#include "TempLat/lattice/algebra/operators/heavisidestepfunction.h"
#include "TempLat/lattice/algebra/operators/log.h"
#include "TempLat/lattice/algebra/operators/multiply.h"
#include "TempLat/lattice/algebra/operators/operators.h"
#include "TempLat/lattice/algebra/operators/power.h"
#include "TempLat/lattice/algebra/operators/shift.h"
#include "TempLat/lattice/algebra/operators/sine.h"
#include "TempLat/lattice/algebra/operators/sinh.h"
#include "TempLat/lattice/algebra/operators/squareroot.h"
#include "TempLat/lattice/algebra/operators/subtract.h"
#include "TempLat/lattice/algebra/operators/tanh.h"
#include "TempLat/lattice/algebra/operators/unaryminus.h"
#include "TempLat/lattice/algebra/operators/unaryoperator.h"

// ---------------------------------------------------------------------------
// Field algebra: operators on lists of fields
// ---------------------------------------------------------------------------
#include "TempLat/lattice/algebra/listoperators/binaryfold.h"
#include "TempLat/lattice/algebra/listoperators/derivatives.h"
#include "TempLat/lattice/algebra/listoperators/foldmultiply.h"
#include "TempLat/lattice/algebra/listoperators/listabsolutevalue.h"
#include "TempLat/lattice/algebra/listoperators/listadd.h"
#include "TempLat/lattice/algebra/listoperators/listbinaryoperator.h"
#include "TempLat/lattice/algebra/listoperators/listcomplexconjugate.h"
#include "TempLat/lattice/algebra/listoperators/listdivide.h"
#include "TempLat/lattice/algebra/listoperators/listexponential.h"
#include "TempLat/lattice/algebra/listoperators/listlaplacian.h"
#include "TempLat/lattice/algebra/listoperators/listlog.h"
#include "TempLat/lattice/algebra/listoperators/listmultiply.h"
#include "TempLat/lattice/algebra/listoperators/listoperators.h"
#include "TempLat/lattice/algebra/listoperators/listpower.h"
#include "TempLat/lattice/algebra/listoperators/listshift.h"
#include "TempLat/lattice/algebra/listoperators/listsquareroot.h"
#include "TempLat/lattice/algebra/listoperators/listsubtract.h"
#include "TempLat/lattice/algebra/listoperators/listunaryminus.h"
#include "TempLat/lattice/algebra/listoperators/listunaryoperator.h"
#include "TempLat/lattice/algebra/listoperators/norm.h"
#include "TempLat/lattice/algebra/listoperators/total.h"
#include "TempLat/lattice/algebra/listoperators/vectordotter.h"

// ---------------------------------------------------------------------------
// Field algebra: conditional getters
// ---------------------------------------------------------------------------
#include "TempLat/lattice/algebra/conditional/conditionalbinarygetter.h"
#include "TempLat/lattice/algebra/conditional/conditionallistbinarygetter.h"
#include "TempLat/lattice/algebra/conditional/conditionallistunarygetter.h"
#include "TempLat/lattice/algebra/conditional/conditionalunarygetter.h"

// ---------------------------------------------------------------------------
// Field algebra: numeric constants and symbols
// ---------------------------------------------------------------------------
#include "TempLat/lattice/algebra/constants/halftype.h"
#include "TempLat/lattice/algebra/constants/number.h"
#include "TempLat/lattice/algebra/constants/numbercollection.h"
#include "TempLat/lattice/algebra/constants/onetype.h"
#include "TempLat/lattice/algebra/constants/symbols.h"
#include "TempLat/lattice/algebra/constants/zerotype.h"

// ---------------------------------------------------------------------------
// Field algebra: coordinates and wavenumbers
// ---------------------------------------------------------------------------
#include "TempLat/lattice/algebra/coordinates/dimensioncountrecorder.h"
#include "TempLat/lattice/algebra/coordinates/momentummultiplicity.h"
#include "TempLat/lattice/algebra/coordinates/spatialcoordinate.h"
#include "TempLat/lattice/algebra/coordinates/wavenumber.h"

// ---------------------------------------------------------------------------
// Field algebra: spatial derivatives
// ---------------------------------------------------------------------------
#include "TempLat/lattice/algebra/spatialderivatives/backdiff.h"
#include "TempLat/lattice/algebra/spatialderivatives/forwdiff.h"
#include "TempLat/lattice/algebra/spatialderivatives/forwdij.h"
#include "TempLat/lattice/algebra/spatialderivatives/latticeforwardgradient.h"
#include "TempLat/lattice/algebra/spatialderivatives/latticelaplacian.h"
#include "TempLat/lattice/algebra/spatialderivatives/neutdiff.h"
#include "TempLat/lattice/algebra/spatialderivatives/neutdij.h"
#include "TempLat/lattice/algebra/spatialderivatives/normgradientsquare.h"

// ---------------------------------------------------------------------------
// Field algebra: complex fields
// ---------------------------------------------------------------------------
#include "TempLat/lattice/algebra/complexalgebra/arg.h"
#include "TempLat/lattice/algebra/complexalgebra/ascomplexfield.h"
#include "TempLat/lattice/algebra/complexalgebra/asfourier.h"
#include "TempLat/lattice/algebra/complexalgebra/complexalgebra.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfield.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfieldadd.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfieldaverager.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfieldbinaryoperator.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfieldconjugate.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfieldfourierview.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfieldmultiply.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfieldoperator.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfieldshift.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfieldsubtract.h"
#include "TempLat/lattice/algebra/complexalgebra/complexfieldunaryoperator.h"
#include "TempLat/lattice/algebra/complexalgebra/complexwrapper.h"
#include "TempLat/lattice/algebra/complexalgebra/helpers/complexfieldget.h"
#include "TempLat/lattice/algebra/complexalgebra/helpers/complexgetgetreturntype.h"
#include "TempLat/lattice/algebra/complexalgebra/helpers/hascomplexfieldget.h"
#include "TempLat/lattice/algebra/complexalgebra/imag.h"
#include "TempLat/lattice/algebra/complexalgebra/real.h"
#include "TempLat/lattice/algebra/complexalgebra/scalarcomplexmultiply.h"

// ---------------------------------------------------------------------------
// Field algebra: SU(2) fields and doublets
// ---------------------------------------------------------------------------
#include "TempLat/lattice/algebra/su2algebra/complexfieldsu2doubletmultiply.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/hassu2doubletget.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/hassu2get.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/paulivectorsalgebra.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/su2doubletget.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/su2doubletgetgetreturntype.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/su2get.h"
#include "TempLat/lattice/algebra/su2algebra/helpers/su2getgetreturntype.h"
#include "TempLat/lattice/algebra/su2algebra/scalarsu2multiplication.h"
#include "TempLat/lattice/algebra/su2algebra/su2algebra.h"
#include "TempLat/lattice/algebra/su2algebra/su2averager.h"
#include "TempLat/lattice/algebra/su2algebra/su2binaryoperator.h"
#include "TempLat/lattice/algebra/su2algebra/su2commutator.h"
#include "TempLat/lattice/algebra/su2algebra/su2dagger.h"
#include "TempLat/lattice/algebra/su2algebra/su2dotter.h"
#include "TempLat/lattice/algebra/su2algebra/su2doublet.h"
#include "TempLat/lattice/algebra/su2algebra/su2doubletaverager.h"
#include "TempLat/lattice/algebra/su2algebra/su2doubletbinaryoperator.h"
#include "TempLat/lattice/algebra/su2algebra/su2doubletdagger.h"
#include "TempLat/lattice/algebra/su2algebra/su2doubletdotter.h"
#include "TempLat/lattice/algebra/su2algebra/su2doubletoperator.h"
#include "TempLat/lattice/algebra/su2algebra/su2doubletshift.h"
#include "TempLat/lattice/algebra/su2algebra/su2doubletsubtract.h"
#include "TempLat/lattice/algebra/su2algebra/su2doubletsum.h"
#include "TempLat/lattice/algebra/su2algebra/su2doubletunaryoperator.h"
#include "TempLat/lattice/algebra/su2algebra/su2doubletwrapper.h"
#include "TempLat/lattice/algebra/su2algebra/su2expmap.h"
#include "TempLat/lattice/algebra/su2algebra/su2expmapinv.h"
#include "TempLat/lattice/algebra/su2algebra/su2field.h"
#include "TempLat/lattice/algebra/su2algebra/su2generators.h"
#include "TempLat/lattice/algebra/su2algebra/su2groupwrapper.h"
#include "TempLat/lattice/algebra/su2algebra/su2liealgebrafield.h"
#include "TempLat/lattice/algebra/su2algebra/su2multiply.h"
#include "TempLat/lattice/algebra/su2algebra/su2operator.h"
#include "TempLat/lattice/algebra/su2algebra/su2shift.h"
#include "TempLat/lattice/algebra/su2algebra/su2su2doubletmultiply.h"
#include "TempLat/lattice/algebra/su2algebra/su2subtract.h"
#include "TempLat/lattice/algebra/su2algebra/su2sum.h"
#include "TempLat/lattice/algebra/su2algebra/su2trace.h"
#include "TempLat/lattice/algebra/su2algebra/su2tracedeficit.h"
#include "TempLat/lattice/algebra/su2algebra/su2unaryoperator.h"
#include "TempLat/lattice/algebra/su2algebra/su2wrapper.h"

// ---------------------------------------------------------------------------
// Field algebra: 3x3 matrix fields
// ---------------------------------------------------------------------------
#include "TempLat/lattice/algebra/matrix3x3algebra/allmatrixcomponents.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/allmatrixoperator.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hashermget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hasmatrixget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hassymget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hassymtracelessget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/hermget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/matrixget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/matrixgetgetreturntype.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/symget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/helpers/symtracelessget.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/hermsymtracelessmultiply.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/hermsymtracelessmultiplytrace.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/hermwrapper.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/matrix3x3algebra.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/matrixbinaryoperator.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/matrixmatrixmultiply.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/matrixmatrixmultiplytrace.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/matrixwrapper.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/scalarsymtracelessmultiply.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symsymtracelessmultiply.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symsymtracelessmultiplytrace.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessadd.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessbinaryoperator.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessconjugate.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessfield.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessfieldasfourier.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessfieldfourierview.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessfieldshift.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelesssubtract.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelessunaryoperator.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symtracelesswrapper.h"
#include "TempLat/lattice/algebra/matrix3x3algebra/symwrapper.h"

// ---------------------------------------------------------------------------
// Field algebra: gauge fields
// ---------------------------------------------------------------------------
#include "TempLat/lattice/algebra/gaugealgebra/backwardcovariantderivative.h"
#include "TempLat/lattice/algebra/gaugealgebra/centeredcovariantderivative.h"
#include "TempLat/lattice/algebra/gaugealgebra/centeredcovariantderivativeo4.h"
#include "TempLat/lattice/algebra/gaugealgebra/centeredfields.h"
#include "TempLat/lattice/algebra/gaugealgebra/covariantlaplacian.h"
#include "TempLat/lattice/algebra/gaugealgebra/fieldstrength.h"
#include "TempLat/lattice/algebra/gaugealgebra/forwardcovariantderivative.h"
#include "TempLat/lattice/algebra/gaugealgebra/gaugealgebra.h"
#include "TempLat/lattice/algebra/gaugealgebra/magneticfield.h"
#include "TempLat/lattice/algebra/gaugealgebra/nonabelianclover.h"
#include "TempLat/lattice/algebra/gaugealgebra/plaquette.h"
#include "TempLat/lattice/algebra/gaugealgebra/plaquetteback.h"
#include "TempLat/lattice/algebra/gaugealgebra/u1exponential.h"
#include "TempLat/lattice/algebra/gaugealgebra/u1field.h"

// ---------------------------------------------------------------------------
// Field algebra: axionic couplings
// ---------------------------------------------------------------------------
#include "TempLat/lattice/algebra/axionalgebra/electricfield2.h"
#include "TempLat/lattice/algebra/axionalgebra/magneticfield4.h"

// ---------------------------------------------------------------------------
// Field algebra: random fields
// ---------------------------------------------------------------------------
#include "TempLat/lattice/algebra/random/randomgaussianfield.h"

// ---------------------------------------------------------------------------
// Field algebra: adapters
// ---------------------------------------------------------------------------
#include "TempLat/lattice/algebra/adapter/momentuminterpolator.h"

// ---------------------------------------------------------------------------
// Field algebra: traits and helpers used to build expressions
// ---------------------------------------------------------------------------
#include "TempLat/lattice/algebra/helpers/confirmghosts.h"
#include "TempLat/lattice/algebra/helpers/confirmspace.h"
#include "TempLat/lattice/algebra/helpers/doeval.h"
#include "TempLat/lattice/algebra/helpers/getcomponent.h"
#include "TempLat/lattice/algebra/helpers/getderiv.h"
#include "TempLat/lattice/algebra/helpers/getdx.h"
#include "TempLat/lattice/algebra/helpers/getfloattype.h"
#include "TempLat/lattice/algebra/helpers/getgetreturntype.h"
#include "TempLat/lattice/algebra/helpers/getkir.h"
#include "TempLat/lattice/algebra/helpers/getndim.h"
#include "TempLat/lattice/algebra/helpers/getngrid.h"
#include "TempLat/lattice/algebra/helpers/getstring.h"
#include "TempLat/lattice/algebra/helpers/gettoolbox.h"
#include "TempLat/lattice/algebra/helpers/getvectorcomponent.h"
#include "TempLat/lattice/algebra/helpers/getvectorsize.h"
#include "TempLat/lattice/algebra/helpers/ghostshunter.h"
#include "TempLat/lattice/algebra/helpers/hasderivmethod.h"
#include "TempLat/lattice/algebra/helpers/hasdoweneedghosts.h"
#include "TempLat/lattice/algebra/helpers/hasdx.h"
#include "TempLat/lattice/algebra/helpers/haseval.h"
#include "TempLat/lattice/algebra/helpers/hasexplicitcoordinatedependence.h"
#include "TempLat/lattice/algebra/helpers/hasghostmethod.h"
#include "TempLat/lattice/algebra/helpers/haskir.h"
#include "TempLat/lattice/algebra/helpers/hasspaceconfirmationmethods.h"
#include "TempLat/lattice/algebra/helpers/hasstaticgetter.h"
#include "TempLat/lattice/algebra/helpers/hasstringmethod.h"
#include "TempLat/lattice/algebra/helpers/hastoolbox.h"
#include "TempLat/lattice/algebra/helpers/hasvectorgetmethod.h"
#include "TempLat/lattice/algebra/helpers/isarithmetic.h"
#include "TempLat/lattice/algebra/helpers/iscomplextype.h"
#include "TempLat/lattice/algebra/helpers/isscalartype.h"
#include "TempLat/lattice/algebra/helpers/isstdgettable.h"
#include "TempLat/lattice/algebra/helpers/istemplatgettable.h"
#include "TempLat/lattice/algebra/helpers/isvariadicindex.h"
#include "TempLat/lattice/algebra/helpers/postget.h"
#include "TempLat/lattice/algebra/helpers/preget.h"

// ---------------------------------------------------------------------------
// Field algebra: remaining headers
// ---------------------------------------------------------------------------
#include "TempLat/lattice/algebra/algebra.h"
#include "TempLat/lattice/algebra/spacestateinterface.h"

// ---------------------------------------------------------------------------
// Measurements: averagers, maxima and radial projections
// ---------------------------------------------------------------------------
#include "TempLat/lattice/measuringtools/accumulatortype.h"
#include "TempLat/lattice/measuringtools/averager.h"
#include "TempLat/lattice/measuringtools/averagerhelper.h"
#include "TempLat/lattice/measuringtools/maximum.h"
#include "TempLat/lattice/measuringtools/measurements.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/kbins.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialbincomputer.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionrebinner.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionresult.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionsinglebinandvalue.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionsingledatum.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/radialprojectionsinglequantity.h"
#include "TempLat/lattice/measuringtools/projectionhelpers/unbinnedradialprojectionresult.h"
#include "TempLat/lattice/measuringtools/radialprojector.h"
#include "TempLat/lattice/measuringtools/toolwithownmemory.h"
#include "TempLat/lattice/measuringtools/wallaverager.h"

// ---------------------------------------------------------------------------
// Input and output
// ---------------------------------------------------------------------------
#include "TempLat/lattice/IO/HDF5/fileloaderhdf5.h"
#include "TempLat/lattice/IO/HDF5/filesaverhdf5.h"
#include "TempLat/lattice/IO/HDF5/helpers/hdf5attribute.h"
#include "TempLat/lattice/IO/HDF5/helpers/hdf5dataset.h"
#include "TempLat/lattice/IO/HDF5/helpers/hdf5file.h"
#include "TempLat/lattice/IO/HDF5/helpers/hdf5group.h"
#include "TempLat/lattice/IO/HDF5/helpers/hdf5object.h"
#include "TempLat/lattice/IO/HDF5/helpers/hdf5timeseries.h"
#include "TempLat/lattice/IO/HDF5/helpers/hdf5type.h"
#include "TempLat/lattice/IO/fileio.h"

// ---------------------------------------------------------------------------
// Utilities: logging
// ---------------------------------------------------------------------------
#include "TempLat/util/log/colors.h"
#include "TempLat/util/log/log.h"
#include "TempLat/util/log/logmutex.h"
#include "TempLat/util/log/puttostream.h"
#include "TempLat/util/log/saycomplete.h"
#include "TempLat/util/log/streamcacher.h"
#include "TempLat/util/log/strippathfromfilename.h"
#include "TempLat/util/log/timespent.h"
#include "TempLat/util/log/trailingzerochar.h"

// ---------------------------------------------------------------------------
// Utilities: conditional output streams
// ---------------------------------------------------------------------------
#include "TempLat/util/conditionaloutput/conditionalfilestream.h"
#include "TempLat/util/conditionaloutput/conditionalsayshort.h"
#include "TempLat/util/conditionaloutput/conditionalstream.h"
#include "TempLat/util/conditionaloutput/outputstream.h"

// ---------------------------------------------------------------------------
// Utilities: compile-time tags and range iteration
// ---------------------------------------------------------------------------
#include "TempLat/util/rangeiteration/for_in_range.h"
#include "TempLat/util/rangeiteration/make_list_tag.h"
#include "TempLat/util/rangeiteration/make_tuple_tag.h"
#include "TempLat/util/rangeiteration/sum_in_range.h"
#include "TempLat/util/rangeiteration/tag.h"
#include "TempLat/util/rangeiteration/taglist.h"
#include "TempLat/util/rangeiteration/tagliteral.h"

// ---------------------------------------------------------------------------
// Utilities: random numbers
// ---------------------------------------------------------------------------
#include "TempLat/util/random/randomgaussian.h"
#include "TempLat/util/random/randomuniform.h"

// ---------------------------------------------------------------------------
// Utilities: hashing
// ---------------------------------------------------------------------------
#include "TempLat/util/hash/keccakhash.h"
#include "TempLat/util/hash/keccakhashbareclass.h"

// ---------------------------------------------------------------------------
// Utilities: JSON
// ---------------------------------------------------------------------------
#include "TempLat/util/json/filetojson.h"
#include "TempLat/util/json/simplejson.h"

// ---------------------------------------------------------------------------
// Utilities: debugging
// ---------------------------------------------------------------------------
#include "TempLat/util/debug/cdemangle.h"
#include "TempLat/util/debug/mpidebuggerhanger.h"
#include "TempLat/util/debug/poormansprofile.h"
#include "TempLat/util/debug/stacktrace.h"
#include "TempLat/util/debug/stacktraceplainptrs.h"
#include "TempLat/util/debug/stacktraceptrtofileaddr.h"

// ---------------------------------------------------------------------------
// Utilities: general
// ---------------------------------------------------------------------------
#include "TempLat/util/almostequal.h"
#include "TempLat/util/assignabletuple.h"
#include "TempLat/util/benchmark.h"
#include "TempLat/util/concat.h"
#include "TempLat/util/constants.h"
#include "TempLat/util/constexpr_for.h"
#include "TempLat/util/containsspace.h"
#include "TempLat/util/cstyletime.h"
#include "TempLat/util/demangle.h"
#include "TempLat/util/endianness.h"
#include "TempLat/util/exception.h"
#include "TempLat/util/factorize.h"
#include "TempLat/util/filetostring.h"
#include "TempLat/util/flattenstd.h"
#include "TempLat/util/flattentuple.h"
#include "TempLat/util/floattostring.h"
#include "TempLat/util/foreach.h"
#include "TempLat/util/function.h"
#include "TempLat/util/getcpptypename.h"
#include "TempLat/util/iscomposite.h"
#include "TempLat/util/isincontainer.h"
#include "TempLat/util/istuplelike.h"
#include "TempLat/util/latinindiceslist.h"
#include "TempLat/util/loadbalance.h"
#include "TempLat/util/makeflatlist.h"
#include "TempLat/util/makeuniformarray.h"
#include "TempLat/util/namedtmpfile.h"
#include "TempLat/util/ndloop.h"
#include "TempLat/util/numericalintegrator.h"
#include "TempLat/util/parenthesisstripper.h"
#include "TempLat/util/powr.h"
#include "TempLat/util/prettytostring.h"
#include "TempLat/util/shiftedindexsequence.h"
#include "TempLat/util/spline.h"
#include "TempLat/util/static_max.h"
#include "TempLat/util/staticif.h"
#include "TempLat/util/staticwarning.h"
#include "TempLat/util/stdatomictype.h"
#include "TempLat/util/stringify.h"
#include "TempLat/util/stringtrimmer.h"
#include "TempLat/util/templatarray.h"
#include "TempLat/util/templatvector.h"
#include "TempLat/util/timer.h"
#include "TempLat/util/tuple_size.h"
#include "TempLat/util/tuple_tools.h"
#include "TempLat/util/tuplemaker.h"

#endif
