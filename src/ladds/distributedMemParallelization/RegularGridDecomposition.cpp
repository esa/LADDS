/**
 * @file RegularGridDecomposition.cpp
 * @author F. Gratl
 * @date 20.05.22
 */

#include "RegularGridDecomposition.h"

#include <autopas/utils/ArrayMath.h>
#include <autopas/utils/WrapOpenMP.h>

#include "ladds/distributedMemParallelization/ParticleCommunicator.h"

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
  printMPIInfo();
}

int LADDS::RegularGridDecomposition::getRank(const std::array<double, 3> &coordinates) const {
  using autopas::utils::ArrayMath::abs;
  using autopas::utils::ArrayMath::div;
  using autopas::utils::ArrayMath::sub;
  using autopas::utils::ArrayUtils::static_cast_array;

  // we need to translate to 0 avoid problems with negative coordinates
  const auto translatedCoords = sub(coordinates, globalBoxMin);
  const auto localBoxLength = abs(sub(localBoxMax, localBoxMin));
  const auto rankGridCoords = static_cast_array<int>(div(translatedCoords, localBoxLength));
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

std::vector<LADDS::Particle> LADDS::RegularGridDecomposition::communicateParticles(
    std::vector<LADDS::Particle> &leavingParticles, autopas::AutoPas<LADDS::Particle> &autopas) const {
  // Set up the communicator
  const auto &comm = getCommunicator();
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

  const auto &[coordsMax, periods, coords] = getGridInfo();
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
}