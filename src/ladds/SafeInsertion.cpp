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
      delayedInsertion.push_back(nSat);
    }
  }
  return delayedInsertion;
}
