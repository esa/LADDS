/**
 * @file ParticleMigrationHandler.cpp
 * @author P. Gomez
 * @date 2022-07-18
 */
#include "ParticleMigrationHandler.h"

#include "ladds/particle/Particle.h"

std::tuple<LADDS::CollisionFunctor::CollisionCollectionT, LADDS::CollisionFunctor::CollisionCollectionT>
LADDS::ParticleMigrationHandler::collisionDetectionImmigrants(AutoPas_t &autopas,
                                                              std::vector<LADDS::Particle> &incomingParticles,
                                                              double deltaT,
                                                              double maxV,
                                                              double collisionDistanceFactor,
                                                              double minDetectionRadius,
                                                              double CDMcutoffInKM) {
  using autopas::utils::ArrayMath::abs;
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::div;
  using autopas::utils::ArrayMath::mul;
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::sub;

  LADDS::CollisionFunctor collisionFunctor(
      autopas.getCutoff(), deltaT, collisionDistanceFactor, minDetectionRadius, CDMcutoffInKM);
  const std::array<double, 3> maxVVec{maxV, maxV, maxV};
  const auto maxCoveredDistance = mulScalar(maxVVec, deltaT);

  for (auto &immigrant : incomingParticles) {
    // find the bounding box around the path the immigrant took since the last update
    // pos - velocity * time
    const auto posFirstPotentialCollision = sub(immigrant.getPosition(), mulScalar(immigrant.getVelocity(), deltaT));
    // pos + maxVelocity * time
    const auto higherCorner = add(posFirstPotentialCollision, maxCoveredDistance);
    // pos - maxVelocity * time
    const auto lowerCorner = sub(posFirstPotentialCollision, maxCoveredDistance);

    // interact the immigrant with everything in the box
    autopas.forEachInRegion(
        [&](auto &p) { collisionFunctor.AoSFunctor(immigrant, p, false); }, lowerCorner, higherCorner);
  }

  return {collisionFunctor.getCollisions(), collisionFunctor.getEvadedCollisions()};
}
