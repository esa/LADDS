/**
 * @file iteratePairwiseCollisionFunctor.cpp
 * @author F. Gratl
 * @date 28.06.21
 *
 * Contains a explicit template instantiation for the iteratePairwise() method with the collision functor of the main
 * AutoPas class and the particle type used here. This is linked into the main executable to enable the
 * other compilation units to only declare, but not instantiate this template.
 */

#include <autopas/AutoPasImpl.h>

#include "ladds/CollisionFunctor.h"
#include "ladds/Debris.h"

//! @cond Doxygen_Suppress
template bool autopas::AutoPas<Debris>::iteratePairwise(CollisionFunctor *);
//! @endcond
