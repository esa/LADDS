/**
 * @file ConjunctionLogger.cpp
 * @author F. Gratl
 * @date 06.12.21
 */

#include <autopas/utils/Timer.h>
#include <spdlog/async.h>

#include "ConjunctionLogger.h"
namespace LADDS {

ConjunctionLogger::ConjunctionLogger(const std::string &outputSuffix) {
  const auto outputFileName("conjunctions_" + outputSuffix + autopas::utils::Timer::getDateStamp() + ".csv");
  // create and register the actual logger
  auto logger = spdlog::basic_logger_mt<spdlog::async_factory>(_loggerName, outputFileName);
  logger->set_level(spdlog::level::info);
  // set pattern to only log what is given
  logger->set_pattern("%v");
  // log the header
  logger->info("Iteration,P1,P2,SquaredDistance");
}

ConjunctionLogger::~ConjunctionLogger() {
  spdlog::drop(_loggerName);
}

void ConjunctionLogger::log(size_t iteration, const Particle &p1, const Particle &p2, double distance) {
  // compute the squared distance inside a lambda inside the macro so it will not be executed when the logger is
  // compile-disabled
  SPDLOG_LOGGER_INFO(spdlog::get(_loggerName), "{},{},{},{}", iteration, p1.getID(), p2.getID(), distance);
}

void ConjunctionLogger::writeConjunctions(size_t iteration, const CollisionFunctor::CollisionCollectionT &collisions) {
  for (const auto &[p1, p2, distanceSquare, _] : collisions) {
    log(iteration, *p1, *p2, distanceSquare);
  }
}

void ConjunctionLogger::writeEvasions(size_t iteration, const CollisionFunctor::CollisionCollectionT &evasions) {
  for (const auto &[p1, p2, distanceSquare, _] : evasions) {
    log(iteration, *p1, *p2, distanceSquare);
  }
}

}  // namespace LADDS