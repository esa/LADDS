/**
 * @file SafeInsertion.h
 * @author albert
 * @date 27.11.21
 */

#pragma once

#include <autopas/AutoPasDecl.h>

#include <iostream>
#include <vector>

#include "ladds/particle/Particle.h"

namespace SafeInsertion {

/**
 * inserts the particles of newSatellites to autopas that have a safe distance to particles
 * in the simulation and returns particles that are not added
 * @param autopas
 * @param newSatellites
 * @param cutoff
 * @return std::vector<Particle> delayedInsertion
 */
std::vector<Particle> insert(autopas::AutoPas<Particle> &autopas,
                             const std::vector<Particle> &newSatellites,
                             double cutoff);

}  // namespace SafeInsertion