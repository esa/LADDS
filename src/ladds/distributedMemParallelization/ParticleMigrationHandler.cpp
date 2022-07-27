/**
 * @file ParticleMigrationHandler.cpp
 * @author P. Gomez
 * @date 2022-07-18
 */
#include "ParticleMigrationHandler.h"

#include "ladds/distributedMemParallelization/ParticleCommunicator.h"
#include "ladds/particle/Particle.h"

std::vector<LADDS::Particle> LADDS::ParticleMigrationHandler::communicateParticles(
    std::vector<LADDS::Particle> &leavingParticles,
    autopas::AutoPas<LADDS::Particle> &autopas,
    const DomainDecomposition &decomposition) {
  // Set up the communicator
  const auto &comm = decomposition.getCommunicator();
  const auto &localBoxMin = autopas.getBoxMin();
  const auto &localBoxMax = autopas.getBoxMax();
  auto logger = spdlog::get(LADDS_SPD_LOGGER_NAME);

  ParticleCommunicator particleCommunicator;
  std::vector<LADDS::Particle> incomingParticles;

  auto getNeighborRank = [&](const auto &coordsThis, int direction, auto op) {
    auto coordsOther = coordsThis;
    coordsOther[direction] = op(coordsOther[direction], 1);
    int rankOther{};
    autopas::AutoPas_MPI_Cart_rank(comm, coordsOther.data(), &rankOther);
    return rankOther;
  };

  if (const auto *regularGridDecomp = dynamic_cast<const RegularGridDecomposition *>(&decomposition)) {
    const auto &[coordsMax, periods, coords] = regularGridDecomp->getGridInfo();
    enum CommDir : int { X, Y, Z };
    for (CommDir commDir = X; commDir <= Z; commDir = static_cast<CommDir>(commDir + 1)) {
      // trigger both non-blocking sends before doing both blocking receives

      // send left (negative direction)
      if (coords[commDir] != 0) {
        // sort particles that are leaving in the negative direction to the end of leavingParticles
        auto leavingParticlesIter =
            std::partition(leavingParticles.begin(), leavingParticles.end(), [&](const Particle &p) {
              return p.getPosition()[commDir] > localBoxMin[commDir];
            });
        const int rankLeft = getNeighborRank(coords, commDir, std::minus<>());
        particleCommunicator.sendParticles(leavingParticlesIter, leavingParticles.end(), rankLeft, comm);
        // clip sent particles
        leavingParticles.erase(leavingParticlesIter, leavingParticles.end());
      }

      // communication right (positive direction)
      if (coords[commDir] != coordsMax[commDir] - 1) {
        // sort particles that are leaving in the positive direction to the end of leavingParticles
        auto leavingParticlesIter =
            std::partition(leavingParticles.begin(), leavingParticles.end(), [&](const Particle &p) {
              return p.getPosition()[commDir] < localBoxMax[commDir];
            });
        const int rankRight = getNeighborRank(coords, commDir, std::plus<>());
        particleCommunicator.sendParticles(leavingParticlesIter, leavingParticles.end(), rankRight, comm);
        // clip sent particles
        leavingParticles.erase(leavingParticlesIter, leavingParticles.end());

        // receive
        auto incomingParticlesRight = particleCommunicator.receiveParticles(rankRight, comm);
        incomingParticles.insert(incomingParticles.end(), incomingParticlesRight.begin(), incomingParticlesRight.end());
      }

      // receive left (negative direction)
      if (coords[commDir] != 0) {
        const int rankLeft = getNeighborRank(coords, commDir, std::minus<>());
        auto incomingParticlesLeft = particleCommunicator.receiveParticles(rankLeft, comm);
        incomingParticles.insert(incomingParticles.end(), incomingParticlesLeft.begin(), incomingParticlesLeft.end());
      }

      // sort particles which have to be added to the local container to the front of incomingParticles
      auto incomingParticlesIter =
          std::partition(incomingParticles.begin(), incomingParticles.end(), [&](const Particle &p) {
            return autopas::utils::inBox(p.getPosition(), localBoxMin, localBoxMax);
          });

      // move particles that were not inserted to leaving particles.
      // These pass-through particles are those changing ranks in more than one dimension.
      leavingParticles.insert(leavingParticles.end(), incomingParticlesIter, incomingParticles.end());
      incomingParticles.erase(incomingParticlesIter, incomingParticles.end());

      particleCommunicator.waitAndFlushBuffers();
    }
    return incomingParticles;
  } else if (const auto *altitudeBasedDecomposition =
                 dynamic_cast<const AltitudeBasedDecomposition *>(&decomposition)) {
    const auto rank = altitudeBasedDecomposition->getRank();
    int numRanks{};
    autopas::AutoPas_MPI_Comm_size(decomposition.getCommunicator(), &numRanks);
    const auto coords = std::array<int, 1>{rank};
    // trigger both non-blocking sends before doing both blocking receives
    // send left (negative direction), commDir is only x for this decomp, thus 0
    const int commDir = 0;
    const auto altBoxMinSquared = std::pow(altitudeBasedDecomposition->getAltitudeOfRank(rank), 2.);
    const auto altBoxMaxSquared = std::pow(altitudeBasedDecomposition->getAltitudeOfRank(rank + 1), 2.);

    if (rank != 0) {
      // sort particles that are leaving in the negative direction to the end of leavingParticles
      auto leavingParticlesIter =
          std::partition(leavingParticles.begin(), leavingParticles.end(), [&](const Particle &p) {
            return autopas::utils::ArrayMath::dot(p.getPosition(), p.getPosition()) > altBoxMinSquared;
          });
      const int rankLeft = getNeighborRank(coords, commDir, std::minus<>());

      // get left neighbors lower limit to check particle is not skipping right through it
      const auto leftNeighboraltBoxMinSquared = std::pow(altitudeBasedDecomposition->getAltitudeOfRank(rankLeft), 2);
      // sanity check
      for (const auto &p : leavingParticles) {
        if (autopas::utils::ArrayMath::dot(p.getPosition(), p.getPosition()) < leftNeighboraltBoxMinSquared) {
          SPDLOG_LOGGER_WARN(
              logger.get(), "Particle {} skipping through left neighbor {}  while migrating.", p.getID(), rankLeft);
        }
      }
      particleCommunicator.sendParticles(leavingParticlesIter, leavingParticles.end(), rankLeft, comm);

      // clip sent particles
      leavingParticles.erase(leavingParticlesIter, leavingParticles.end());
    }

    // communication right (positive direction)
    if (rank != numRanks - 1) {
      // sort particles that are leaving in the positive direction to the end of leavingParticles
      auto leavingParticlesIter =
          std::partition(leavingParticles.begin(), leavingParticles.end(), [&](const Particle &p) {
            return autopas::utils::ArrayMath::dot(p.getPosition(), p.getPosition()) > altBoxMaxSquared;
          });
      const int rankRight = getNeighborRank(coords, commDir, std::plus<>());

      // get right neighbors upper limit to check particle is not skipping right through it
      const auto leftNeighboraltBoxMaxSquared = std::pow(altitudeBasedDecomposition->getAltitudeOfRank(rankRight), 2);
      for (auto &p : leavingParticles) {
        if (autopas::utils::ArrayMath::dot(p.getPosition(), p.getPosition()) > leftNeighboraltBoxMaxSquared) {
          SPDLOG_LOGGER_WARN(
              logger.get(), "Particle {} skipping through right neighbor {}  while migrating.", p.getID(), rankRight);
        }
      }

      particleCommunicator.sendParticles(leavingParticlesIter, leavingParticles.end(), rankRight, comm);

      // clip sent particles
      leavingParticles.erase(leavingParticlesIter, leavingParticles.end());

      // receive
      auto incomingParticlesRight = particleCommunicator.receiveParticles(rankRight, comm);
      incomingParticles.insert(incomingParticles.end(), incomingParticlesRight.begin(), incomingParticlesRight.end());
    }

    // receive left (negative direction)
    if (rank != 0) {
      const int rankLeft = getNeighborRank(coords, commDir, std::minus<>());
      auto incomingParticlesLeft = particleCommunicator.receiveParticles(rankLeft, comm);
      incomingParticles.insert(incomingParticles.end(), incomingParticlesLeft.begin(), incomingParticlesLeft.end());
    }

    leavingParticles.erase(leavingParticles.begin(), leavingParticles.end());
    particleCommunicator.waitAndFlushBuffers();

    SPDLOG_LOGGER_DEBUG(logger.get(), "Rank {} received {} particles", rank, incomingParticles.size());
    return incomingParticles;
  } else {
    throw std::runtime_error("No particle communication implemented for the chosen decomposition!");
  }
}

LADDS::CollisionFunctor::CollisionCollectionT LADDS::ParticleMigrationHandler::collisionDetectionImmigrants(
    AutoPas_t &autopas,
    std::vector<LADDS::Particle> &incomingParticles,
    double deltaT,
    double maxV,
    double collisionDistanceFactor,
    double minDetectionRadius) {
  using autopas::utils::ArrayMath::abs;
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::div;
  using autopas::utils::ArrayMath::mul;
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::sub;

  LADDS::CollisionFunctor collisionFunctor(autopas.getCutoff(), deltaT, collisionDistanceFactor, minDetectionRadius);
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

  return collisionFunctor.getCollisions();
}
