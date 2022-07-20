/**
 * @file AltitudeBasedDecomposition.cpp
 * @author P. Gomez
 * @date 14.07.22
 */

#include "AltitudeBasedDecomposition.h"

#include <autopas/utils/ArrayMath.h>
#include <autopas/utils/WrapOpenMP.h>

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
  // TODO: find a more even distribution of the particles
  altitudeIntervals = this->logspace(std::log10(minAltitude), std::log10(maxAltitude - Physics::R_EARTH), numRanks);
  for (int i = 0; i < numRanks; i++) {
    altitudeIntervals[i] += Physics::R_EARTH;
  }
  altitudeIntervals[0] = 0.0;
  // note that actual max altitude is the corner of the box, so higher than maxaltitude in the config but there should
  // be very satellites between maxAltitude and this
  altitudeIntervals.push_back(std::sqrt(3 * std::pow(maxAltitude, 2.0)));
  SPDLOG_LOGGER_INFO(config.getLogger().get(),
                     "Computed altitude intervals fo ranks: {}",
                     autopas::utils::ArrayUtils::to_string(altitudeIntervals));

  localBoxMin = {-altitudeIntervals[rank + 1], -altitudeIntervals[rank + 1], -altitudeIntervals[rank + 1]};
  localBoxMax = {altitudeIntervals[rank + 1], altitudeIntervals[rank + 1], altitudeIntervals[rank + 1]};

  // print parallelization info
  if (rank == 0) {
    SPDLOG_LOGGER_INFO(config.getLogger().get(),
                       "Parallelization Configuration\n"
                       "MPI Ranks              : {}\n"
                       "OpenMP Threads per Rank: {}\n"
                       "MPI Decomposition      : {}\n",
                       numRanks,
                       autopas::autopas_get_max_threads(),
                       autopas::utils::ArrayUtils::to_string(dims));
  }
}

int LADDS::AltitudeBasedDecomposition::getRank(const std::array<double, 3> &coordinates) const {
  // std::cout << "Getting rank for coordinates " << autopas::utils::ArrayUtils::to_string(coordinates) << ::std::endl;

  using autopas::utils::ArrayMath::div;
  using autopas::utils::ArrayMath::sub;

  const auto altitudeSquared = autopas::utils::ArrayMath::dot(coordinates, coordinates);
  for (size_t i = 0; i < altitudeIntervals.size() - 1; i++) {
    if (altitudeSquared >= altitudeIntervals[i] * altitudeIntervals[i] &&
        altitudeSquared < altitudeIntervals[i + 1] * altitudeIntervals[i + 1]) {
      return i;
    }
  }
  throw std::runtime_error("Could not find rank for coordinates " + autopas::utils::ArrayUtils::to_string(coordinates));
}

std::tuple<std::array<int, 1>, std::array<int, 1>, std::array<int, 1>> LADDS::AltitudeBasedDecomposition::getGridInfo()
    const {
  std::array<int, 1> dims{};
  std::array<int, 1> periods{};
  std::array<int, 1> coords{};
  autopas::AutoPas_MPI_Cart_get(communicator, dims.size(), dims.data(), periods.data(), coords.data());
  return {dims, periods, coords};
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

std::vector<LADDS::Particle> LADDS::AltitudeBasedDecomposition::getLeavingParticles(const AutoPas_t &autopas) const {
  std::vector<LADDS::Particle> particles;
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(communicator, &rank);
  for (auto &particle : autopas) {
    if (this->getRank(particle.getPosition()) != rank) {
      particles.push_back(particle);
    }
  }

  // std::cout << "Leaving particles: " << particles.size() << std::endl;
  // for (auto &particle : particles) {
  //   std::cout << "Leaving particle: " << particle.getID() << " "
  //             << autopas::utils::ArrayUtils::to_string(particle.getPosition()) << std::endl;
  // }

  return particles;
}