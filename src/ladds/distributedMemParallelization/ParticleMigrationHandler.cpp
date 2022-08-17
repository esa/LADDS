/**
 * @file ParticleMigrationHandler.cpp
 * @author P. Gomez
 * @date 2022-07-18
 */
#include "ParticleMigrationHandler.h"

#include "ladds/particle/Particle.h"

std::tuple<LADDS::CollisionFunctor::CollisionCollectionT, LADDS::CollisionFunctor::CollisionCollectionT>
LADDS::ParticleMigrationHandler::collisionDetectionAroundParticles(AutoPas_t &autopas,
                                                                   std::vector<LADDS::Particle> &particles,
                                                                   double deltaT,
                                                                   double maxV,
                                                                   double collisionDistanceFactor,
                                                                   double minDetectionRadius,
                                                                   double CDMcutoffInKM,
                                                                   bool checkForInternalCollisions) {
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

  collisionFunctor.initTraversal();
  for (auto &particle : particles) {
    // find the bounding box around the path the particle took since the last update
    // pos - velocity * time
    const auto posFirstPotentialCollision = sub(particle.getPosition(), mulScalar(particle.getVelocity(), deltaT));
    // pos + maxVelocity * time
    const auto higherCorner = add(posFirstPotentialCollision, maxCoveredDistance);
    // pos - maxVelocity * time
    const auto lowerCorner = sub(posFirstPotentialCollision, maxCoveredDistance);

    // interact the particle with everything in the box
    autopas.forEachInRegion(
        [&](auto &p) { collisionFunctor.AoSFunctor(particle, p, false); }, lowerCorner, higherCorner);

    // interact the particle with all other particles in this vector if desired
    if (checkForInternalCollisions) {
      for (auto &particle2 : particles) {
        // Greater than to make sure that the particle is not interacting with itself and once with the other
        if (particle.getID() > particle2.getID()) collisionFunctor.AoSFunctor(particle, particle2, false);
      }
    }
  }
  collisionFunctor.endTraversal(false);

  return {collisionFunctor.getCollisions(), collisionFunctor.getEvadedCollisions()};
}
