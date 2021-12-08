/**
 * @file Timers.h
 * @author F. Gratl
 * @date 06.12.21
 */

#pragma once

#include <autopas/utils/Timer.h>
#include <yaml-cpp/yaml.h>

class Timers {
 public:
  autopas::utils::Timer total{};
  autopas::utils::Timer initialization{};
  autopas::utils::Timer simulation{};
  autopas::utils::Timer integrator{};
  autopas::utils::Timer constellationInsertion{};
  autopas::utils::Timer collisionDetection{};
  autopas::utils::Timer containerUpdate{};
  autopas::utils::Timer output{};

  /**
   * Print timer information to stdout.
   * @param config
   */
  void printTimers(const YAML::Node &config) const;

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
   * Floating point precision for command line output.
   */
  static constexpr int floatStringPrecision = 3;
};