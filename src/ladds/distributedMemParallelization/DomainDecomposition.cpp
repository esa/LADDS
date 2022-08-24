/**
 * @file DomainDecomposition.cpp
 * @author F. Gratl
 * @date 24.05.22
 */
#include "DomainDecomposition.h"

#include <autopas/utils/WrapOpenMP.h>
#include <autopas/utils/inBox.h>
#include <spdlog/spdlog.h>

std::array<double, 3> LADDS::DomainDecomposition::getGlobalBoxMin() const {
  return globalBoxMin;
}
std::array<double, 3> LADDS::DomainDecomposition::getGlobalBoxMax() const {
  return globalBoxMax;
}
std::array<double, 3> LADDS::DomainDecomposition::getLocalBoxMin() const {
  return localBoxMin;
}
std::array<double, 3> LADDS::DomainDecomposition::getLocalBoxMax() const {
  return localBoxMax;
}
bool LADDS::DomainDecomposition::containsGlobal(const std::array<double, 3> &coordinates) const {
  return autopas::utils::inBox(coordinates, globalBoxMin, globalBoxMax);
}
bool LADDS::DomainDecomposition::containsLocal(const std::array<double, 3> &coordinates) const {
  return autopas::utils::inBox(coordinates, localBoxMin, localBoxMax);
}
autopas::AutoPas_MPI_Comm LADDS::DomainDecomposition::getCommunicator() const {
  return communicator;
}

void LADDS::DomainDecomposition::printMPIInfo() const {
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(this->getCommunicator(), &rank);
  if (rank == 0) {
    auto logger = spdlog::get(LADDS_SPD_LOGGER_NAME);
    int numRanks{};
    autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numRanks);
    SPDLOG_LOGGER_INFO(logger.get(),
                       "Parallelization Configuration\n"
                       "MPI Ranks              : {}\n"
                       "OpenMP Threads per Rank: {}\n",
                       numRanks,
                       autopas::autopas_get_max_threads());
  }
}