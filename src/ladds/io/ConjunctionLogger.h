/**
 * @file ConjunctionLogger.h
 * @author F. Gratl
 * @date 06.12.21
 */

#pragma once

#include <string>
#include "ladds/particle/Particle.h"
/**
 * Helper to log data for individual conjunctions to a csv file for easier analysis.
 *
 * It uses an asynchronous spd logger to write a csv file named "conjunctions_<dateStamp>.csv".
 */
class ConjunctionLogger {
 public:
  /**
   * Constructor initializes the logger and sets the output file name.
   * @param outputSuffix Suffix for all output files produced by this class.
   */
  explicit ConjunctionLogger(const std::string &outputSuffix = "");

  /**
   * Destructor drops the logger from the spd registry.
   */
  ~ConjunctionLogger();

  /**
   * Log the given arguments and the internal buffer to the csv file.
   * @param p1
   * @param p2
   */
  void log(size_t iteration, const Particle &p1, const Particle &p2);

 private:
  std::string _loggerName{"conjunctionLogger"};
};
