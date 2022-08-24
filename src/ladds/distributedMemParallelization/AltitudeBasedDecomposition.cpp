/**
 * @file AltitudeBasedDecomposition.cpp
 * @author P. Gomez
 * @date 14.07.22
 */

#include "AltitudeBasedDecomposition.h"

#include <autopas/utils/ArrayMath.h>
#include <autopas/utils/WrapOpenMP.h>
#include <spdlog/spdlog.h>

#include "ladds/TypeDefinitions.h"
#include "ladds/distributedMemParallelization/ParticleCommunicator.h"
#include "satellitePropagator/physics/Constants.h"

LADDS::AltitudeBasedDecomposition::AltitudeBasedDecomposition(LADDS::ConfigReader &config) {
  // initialize a one-dimensional MPI world in which split the altitudes
  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numRanks);
  // let MPI decide a dimension splitting with no periodicity
  std::array<int, 1> periods{false};
  std::array<int, 1> dims{};
  autopas::AutoPas_MPI_Dims_create(numRanks, dims.size(), dims.data());
  // create a communicator
  autopas::AutoPas_MPI_Cart_create(
      AUTOPAS_MPI_COMM_WORLD, dims.size(), dims.data(), periods.data(), true, &communicator);
  // get coordinates of this rank within the decomposition grid
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(communicator, &rank);
  std::array<int, 1> coords{};
  autopas::AutoPas_MPI_Cart_coords(communicator, rank, dims.size(), coords.data());

  // initialize simulation specific data
  const auto maxAltitude = config.get<double>("sim/maxAltitude");
  const auto minAltitude = config.get<double>("sim/minAltitude");

  globalBoxMin = {-maxAltitude, -maxAltitude, -maxAltitude};
  globalBoxMax = {maxAltitude, maxAltitude, maxAltitude};

  // calculate local box extent
  // we have to account for the fact that there won't be anything between 0, Physics::R_EARTH  + minAltitude, so
  // shifting the interval a bit
  altitudeIntervals = this->logspace(std::log10(minAltitude), std::log10(maxAltitude - Physics::R_EARTH), numRanks);
  for (auto &interval : altitudeIntervals) {
    interval += Physics::R_EARTH;
  }
  altitudeIntervals[0] = 0.0;
  // note that actual max altitude is the corner of the box, so higher than maxaltitude in the config but there should
  // be very few satellites between maxAltitude and this
  // maximumAltitude = sqrt(maxAltitude**2 + maxAltitude**2 + maxAltitude**2)
  altitudeIntervals.push_back(std::sqrt(3. * std::pow(maxAltitude, 2.0)));
  SPDLOG_LOGGER_INFO(config.getLogger().get(),
                     "Computed altitude intervals fo ranks: {}",
                     autopas::utils::ArrayUtils::to_string(altitudeIntervals));

  localBoxMin = {-altitudeIntervals[rank + 1], -altitudeIntervals[rank + 1], -altitudeIntervals[rank + 1]};
  localBoxMax = {altitudeIntervals[rank + 1], altitudeIntervals[rank + 1], altitudeIntervals[rank + 1]};

  // print parallelization info
  printMPIInfo();
}

int LADDS::AltitudeBasedDecomposition::getRank(const std::array<double, 3> &coordinates) const {
  auto logger = spdlog::get(LADDS_SPD_LOGGER_NAME);

  using autopas::utils::ArrayMath::div;
  using autopas::utils::ArrayMath::sub;

  const auto altitudeSquared = autopas::utils::ArrayMath::dot(coordinates, coordinates);
  for (size_t i = 0; i < altitudeIntervals.size() - 1; i++) {
    if (altitudeSquared >= altitudeIntervals[i] * altitudeIntervals[i] and
        altitudeSquared < altitudeIntervals[i + 1] * altitudeIntervals[i + 1]) {
      SPDLOG_LOGGER_TRACE(logger.get(),
                          "Getting rank for coordinates {} result was {}",
                          autopas::utils::ArrayUtils::to_string(coordinates),
                          i);
      return i;
    }
  }

  throw std::runtime_error("Could not find rank for coordinates " + autopas::utils::ArrayUtils::to_string(coordinates));
}

int LADDS::AltitudeBasedDecomposition::getRank() const {
  std::array<int, 1> _{};
  std::array<int, 1> __{};
  std::array<int, 1> coords{};
  autopas::AutoPas_MPI_Cart_get(communicator, _.size(), _.data(), __.data(), coords.data());
  return coords[0];
}
std::vector<double> LADDS::AltitudeBasedDecomposition::logspace(const double a, const double b, const int k) {
  double realStart = pow(10., a);
  double realBase = pow(10., (b - a) / k);

  std::vector<double> retval;
  retval.reserve(k);
  for (int i = 0; i < k; i++) {
    retval.push_back(realStart);
    realStart *= realBase;
  }

  return retval;
}

double LADDS::AltitudeBasedDecomposition::getAltitudeOfRank(const int rank) const {
  return altitudeIntervals[rank + 1];
}

std::vector<LADDS::Particle> LADDS::AltitudeBasedDecomposition::getAndRemoveLeavingParticles(AutoPas_t &autopas) const {
  std::vector<LADDS::Particle> particles;
  particles.reserve(autopas.getNumberOfParticles());
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(communicator, &rank);
  auto logger = spdlog::get(LADDS_SPD_LOGGER_NAME);
  for (auto &particle : autopas) {
    SPDLOG_LOGGER_TRACE(logger.get(),
                        "Checking particle {} at Position {} - supposed to be in rank {}",
                        particle.getID(),
                        autopas::utils::ArrayUtils::to_string(particle.getPosition),
                        getRank(particle.getPosition()));
    if (this->getRank(particle.getPosition()) != rank) {
      particles.push_back(particle);
      autopas.deleteParticle(particle);
    }
  }

  return particles;
}

void LADDS::AltitudeBasedDecomposition::rebalanceDecomposition(const std::vector<LADDS::Particle> &particles,
                                                               AutoPas_t &autopas) {
  int numRanks{};
  auto logger = spdlog::get(LADDS_SPD_LOGGER_NAME);
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numRanks);
  if (particles.size() < (unsigned)numRanks) {
    SPDLOG_LOGGER_DEBUG(
        logger.get(),
        "Skipped planned decomposition rebelancing as particle number was lower than amount of ranks: {} vs {}",
        particles.size(),
        numRanks);
    return;
  }
  std::vector<double> squaredAltitudes;
  squaredAltitudes.reserve(particles.size());

  int numberOfBins{};
  int numberOfParticles = particles.size();
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numberOfBins);
  int particlesPerRank = std::ceil(numberOfParticles / numberOfBins);

  // Collect all altitudes
  for (auto &particle : particles) {
    squaredAltitudes.push_back(autopas::utils::ArrayMath::dot(particle.getPosition(), particle.getPosition()));
  }

  // Sort the array, then we can just get the equal partitions based on altitude from that
  std::sort(squaredAltitudes.begin(), squaredAltitudes.end());

  // Now store the altitudes as boundaries, first entry will always be 0, last always remains maxAltitude
  for (int i = 1; i < numberOfBins; i++) {
    // include particle at that position as well
    double altitude = std::sqrt(squaredAltitudes[(i * particlesPerRank) - 1]);

    altitudeIntervals[i] = altitude;
  }

  // get coordinates of this rank within the decomposition grid
  const auto communicator = this->getCommunicator();
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(communicator, &rank);
  localBoxMin = {-altitudeIntervals[rank + 1], -altitudeIntervals[rank + 1], -altitudeIntervals[rank + 1]};
  localBoxMax = {altitudeIntervals[rank + 1], altitudeIntervals[rank + 1], altitudeIntervals[rank + 1]};

  autopas.resizeBox(localBoxMin, localBoxMax);

  SPDLOG_LOGGER_DEBUG(
      logger.get(), "Recomputed altitude intervals: {}", autopas::utils::ArrayUtils::to_string(altitudeIntervals));
}

std::vector<LADDS::Particle> LADDS::AltitudeBasedDecomposition::communicateParticles(
    std::vector<LADDS::Particle> &leavingParticles, autopas::AutoPas<LADDS::Particle> &autopas) const {
  // Set up the communicator
  const auto &comm = getCommunicator();
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

  const auto rank = getRank();
  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(getCommunicator(), &numRanks);
  const auto coords = std::array<int, 1>{rank};
  // trigger both non-blocking sends before doing both blocking receives
  // send left (negative direction), commDir is only x for this decomp, thus 0
  const int commDir = 0;
  const auto altBoxMinSquared = std::pow(getAltitudeOfRank(rank), 2.);
  const auto altBoxMaxSquared = std::pow(getAltitudeOfRank(rank + 1), 2.);

  if (rank != 0) {
    // sort particles that are leaving in the negative direction to the end of leavingParticles
    auto leavingParticlesIter =
        std::partition(leavingParticles.begin(), leavingParticles.end(), [&](const Particle &p) {
          return autopas::utils::ArrayMath::dot(p.getPosition(), p.getPosition()) > altBoxMinSquared;
        });
    const int rankLeft = getNeighborRank(coords, commDir, std::minus<>());

    // get left neighbors lower limit to check particle is not skipping right through it
    const auto leftNeighboraltBoxMinSquared = std::pow(getAltitudeOfRank(rankLeft), 2);
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
    const auto leftNeighboraltBoxMaxSquared = std::pow(getAltitudeOfRank(rankRight), 2);
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
}