/**
 * @file Logger.cpp
 * @author F. Gratl
 * @date 28.06.21
 */

#include "Logger.h"

#include <autopas/utils/WrapMPI.h>

namespace LADDS {

Logger::Logger(std::string name, std::ostream &ostream) : _name(std::move(name)) {
  auto sink = std::make_shared<spdlog::sinks::ostream_sink_mt>(ostream);
  auto logger = std::make_shared<spdlog::logger>(_name, sink);
  // set all loggers which are not on rank 0 to only output errors!
  int mpiRank{};
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &mpiRank);
  int mpiNumRanks{};
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &mpiNumRanks);
  if (mpiRank != 0) {
    logger->set_level(spdlog::level::err);
  }
  // rebuild default format and inject rank number
  std::string mpiRankString = std::to_string(mpiRank);
  std::string mpiNumRanksString = std::to_string(mpiNumRanks);
  mpiRankString.insert(0, mpiNumRanksString.length() - mpiRankString.length(), ' ');
  logger->set_pattern("[%Y-%m-%d %T.%e] [%n] [Rank " + mpiRankString + "] [%l] %v");
  spdlog::register_logger(logger);
}

Logger::~Logger() {
  spdlog::drop(_name);
}

std::shared_ptr<spdlog::logger> Logger::get() const {
  return spdlog::get(_name);
}

const std::string &Logger::getName() const {
  return _name;
}

}  // namespace LADDS