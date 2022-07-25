/**
 * @file RegularGridDecomposition.cpp
 * @author F. Gratl
 * @date 20.05.22
 */

#include "RegularGridDecomposition.h"

#include <autopas/utils/ArrayMath.h>
#include <autopas/utils/WrapOpenMP.h>

LADDS::RegularGridDecomposition::RegularGridDecomposition(LADDS::ConfigReader &config) {
  // initialize MPI
  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numRanks);
  // let MPI decide a dimension splitting with no periodicity
  std::array<int, 3> periods{false, false, false};
  std::array<int, 3> dims{};
  autopas::AutoPas_MPI_Dims_create(numRanks, dims.size(), dims.data());
  // create a communicator
  autopas::AutoPas_MPI_Cart_create(
      AUTOPAS_MPI_COMM_WORLD, dims.size(), dims.data(), periods.data(), true, &communicator);
  // get coordinates of this rank within the decomposition grid
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(communicator, &rank);
  std::array<int, 3> coords{};
  autopas::AutoPas_MPI_Cart_coords(communicator, rank, dims.size(), coords.data());

  // initialize simulation specific data
  const auto maxAltitude = config.get<double>("sim/maxAltitude");
  globalBoxMin = {-maxAltitude, -maxAltitude, -maxAltitude};
  globalBoxMax = {maxAltitude, maxAltitude, maxAltitude};
  // using for readability
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::div;
  using autopas::utils::ArrayMath::mul;
  using autopas::utils::ArrayMath::sub;
  using autopas::utils::ArrayUtils::static_cast_array;
  // calculate local box extent
  const auto globalBoxLength = sub(globalBoxMax, globalBoxMin);
  const auto localBoxLength = div(globalBoxLength, static_cast_array<double>(dims));
  localBoxMin = add(globalBoxMin, mul(localBoxLength, static_cast_array<double>(coords)));
  localBoxMax = add(localBoxMin, localBoxLength);

  // print parallelization info
  if (rank == 0) {
    SPDLOG_LOGGER_INFO(config.getLogger().get(),
                       "Parallelization Configuration\n"
                       "MPI Ranks              : {}\n"
                       "OpenMP Threads per Rank: {}\n"
                       "MPI Decomposition      : {}",
                       numRanks,
                       autopas::autopas_get_max_threads(),
                       autopas::utils::ArrayUtils::to_string(dims));
  }
}

int LADDS::RegularGridDecomposition::getRank(const std::array<double, 3> &coordinates) const {
  using autopas::utils::ArrayMath::div;
  using autopas::utils::ArrayMath::sub;
  using autopas::utils::ArrayUtils::static_cast_array;

  const auto localBoxLength = sub(localBoxMax, localBoxMin);
  const auto rankGridCoords = static_cast_array<int>(div(coordinates, localBoxLength));
  int targetRank{};
  autopas::AutoPas_MPI_Cart_rank(communicator, rankGridCoords.data(), &targetRank);
  return targetRank;
}

std::tuple<std::array<int, 3>, std::array<int, 3>, std::array<int, 3>> LADDS::RegularGridDecomposition::getGridInfo()
    const {
  std::array<int, 3> dims{};
  std::array<int, 3> periods{};
  std::array<int, 3> coords{};
  autopas::AutoPas_MPI_Cart_get(communicator, dims.size(), dims.data(), periods.data(), coords.data());
  return {dims, periods, coords};
}