//
// Created by albert on 27.11.21.
//
#pragma once

#include <autopas/AutoPasDecl.h>

#include <iostream>

#include "Particle.h"

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
                             const double collisionRadius) {
  std::vector<Particle> delayedInsertion = {};

  const double collisionRadiusSquared = collisionRadius * collisionRadius;
  const std::array<double, 3> boxSpan = {collisionRadius, collisionRadius, collisionRadius};
  for (const auto &nSat : newSatellites) {
    // only insert satellites, if they have a reasonable distance to other satellites
    bool collisionFree = true;
    const std::array<double, 3> lowCorner = autopas::utils::ArrayMath::sub(nSat.getPosition(), boxSpan);
    const std::array<double, 3> highCorner = autopas::utils::ArrayMath::add(nSat.getPosition(), boxSpan);
    for (auto iter = autopas.getRegionIterator(lowCorner, highCorner); iter.isValid() and collisionFree; ++iter) {
      std::array<double, 3> diff = autopas::utils::ArrayMath::sub(nSat.getPosition(), iter->getPosition());
      collisionFree = autopas::utils::ArrayMath::dot(diff, diff) > collisionRadiusSquared;
    }
    if (collisionFree) {
      autopas.addParticle(nSat);
    } else {
      std::cout << "delayed!" << std::endl;
      delayedInsertion.push_back(nSat);
    }
  }
  return delayedInsertion;
}

}  // namespace SafeInsertion