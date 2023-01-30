/**
 * @file Timers.h
 * @author F. Gratl
 * @date 06.12.21
 */

#pragma once

#include <autopas/utils/Timer.h>
#include <yaml-cpp/yaml.h>

#include "ConfigReader.h"
#include "ladds/distributedMemParallelization/DomainDecomposition.h"

namespace LADDS {

class Timers {
 public:
  autopas::utils::Timer total{};
  autopas::utils::Timer initialization{};
  autopas::utils::Timer simulation{};
  autopas::utils::Timer integrator{};
  autopas::utils::Timer burnUps{};
  autopas::utils::Timer constellationInsertion{};
  autopas::utils::Timer particleCommunication{};
  autopas::utils::Timer collisionDetection{};
  autopas::utils::Timer collisionDetectionImmigrants{};
  autopas::utils::Timer collisionDetectionEmigrants{};
  autopas::utils::Timer collisionWriting{};
  autopas::utils::Timer collisionSimulation{};
  autopas::utils::Timer evasionWriting{};
  autopas::utils::Timer containerUpdate{};
  autopas::utils::Timer output{};

  /**
   * Print timer information to stdout.
   * @param config
   * @param decomp
   */
  void printTimers(ConfigReader &config, const DomainDecomposition &decomp) const;

 private:
  /**
   * Helper function to pretty-print timers.
   * @param name
   * @param timeNS
   * @param numberWidth
   * @param maxTime
   * @return
   */
  static std::string timerToString(const std::string &name, long timeNS, int numberWidth = 0, long maxTime = 0ul);

  /**
   * Calculate some load imbalance statistics for a given timer
   * @param timer
   * @param decomp
   * @return
   */
  static std::tuple<double, double> calcImbalances(const autopas::utils::Timer &timer,
                                                   const DomainDecomposition &decomp);

  /**
   * Floating point precision for command line output.
   */
  static constexpr int floatStringPrecision = 3;
};

}  // namespace LADDS