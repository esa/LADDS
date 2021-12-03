/**
 * @file SafeInsertion.cpp
 * @author Fabio Gratl
 * @date 03.12.21
 */

#include "SafeInsertion.h"

std::vector<Particle> SafeInsertion::insert(autopas::AutoPas<Particle> &autopas,
                                            const std::vector<Particle> &newSatellites,
                                            double cutoff) {
  std::vector<Particle> delayedInsertion = {};

  const double collisionRadius = 2 * cutoff;
  const double collisionRadiusSquared = collisionRadius * collisionRadius;
  const std::array<double, 3> boxSpan = {collisionRadius, collisionRadius, collisionRadius};
  for (const auto &satellite : newSatellites) {
    // only insert satellites, if they have a reasonable distance to other satellites
    bool collisionFree = true;
    const auto lowCorner = autopas::utils::ArrayMath::sub(satellite.getPosition(), boxSpan);
    const auto highCorner = autopas::utils::ArrayMath::add(satellite.getPosition(), boxSpan);
    for (auto iter = autopas.getRegionIterator(lowCorner, highCorner); iter.isValid() and collisionFree; ++iter) {
      const auto diff = autopas::utils::ArrayMath::sub(satellite.getPosition(), iter->getPosition());
      collisionFree = autopas::utils::ArrayMath::dot(diff, diff) > collisionRadiusSquared;
    }
    if (collisionFree) {
      autopas.addParticle(satellite);
    } else {
      delayedInsertion.push_back(satellite);
    }
  }
  return delayedInsertion;
}
