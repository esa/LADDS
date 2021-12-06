/**
 * @file AutoPasClass.cpp
 * @author F. Gratl
 * @date 28.06.21
 *
 * Contains a explicit template instantiation for the main AutoPas class and the particle type used by md-flexible. This
 * is linked into the main executable to enable the other compilation units to only declare, but not instantiate
 * this template.
 */

#include "autopas/AutoPasImpl.h"
#include "ladds/particle/Particle.h"

template class autopas::AutoPas<Particle>;