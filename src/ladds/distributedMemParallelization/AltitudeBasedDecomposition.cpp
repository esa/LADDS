/**
 * @file AltitudeBasedDecomposition.cpp
 * @author P. Gomez
 * @date 14.07.22
 */

#include "AltitudeBasedDecomposition.h"

#include <autopas/utils/ArrayMath.h>
#include <autopas/utils/WrapOpenMP.h>

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
  // have to add 6371 to avoid ranks inside Earth
  altitude_intervals = this->logspace(std::log10(6371.0 + minAltitude), std::log10(maxAltitude), numRanks);
  altitude_intervals[0] = 0.0;
  altitude_intervals.push_back(maxAltitude);
  SPDLOG_LOGGER_DEBUG(config.getLogger().get(),
                      "Computed altitude intervals: {}",
                      autopas::utils::ArrayUtils::to_string(altitude_intervals));

  localBoxMin = {-altitude_intervals[rank + 1], -altitude_intervals[rank + 1], -altitude_intervals[rank + 1]};
  localBoxMax = {altitude_intervals[rank + 1], altitude_intervals[rank + 1], altitude_intervals[rank + 1]};

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
  std::cout << "Getting rank for coordinates " << autopas::utils::ArrayUtils::to_string(coordinates) << ::std::endl;

  using autopas::utils::ArrayMath::div;
  using autopas::utils::ArrayMath::sub;

  const auto altitudeSquared = autopas::utils::ArrayMath::dot(coordinates, coordinates);

  // all boxes are square
  const auto localBoxLength = localBoxMax[0] - localBoxMin[0];
  const auto location = static_cast<int>(altitudeSquared / (localBoxLength * localBoxLength));
  std::array<int, 1> rankGridCoords{location};
  int targetRank{};
  autopas::AutoPas_MPI_Cart_rank(communicator, rankGridCoords.data(), &targetRank);
  return targetRank;
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