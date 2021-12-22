/**
 * @file ConjunctionLogger.h
 * @author F. Gratl
 * @date 06.12.21
 */

#pragma once

#include <string>

#include "ConjunctionWriterInterface.h"
#include "ladds/particle/Particle.h"
/**
 * Helper to log data for individual conjunctions to a csv file for easier analysis.
 *
 * It uses an asynchronous spd logger to write a csv file named "conjunctions_<dateStamp>.csv".
 */
class ConjunctionLogger final : public ConjuctionWriterInterface {
 public:
  /**
   * Constructor initializes the logger and sets the output file name.
   * @param outputSuffix Suffix for all output files produced by this class.
   */
  explicit ConjunctionLogger(const std::string &outputSuffix = "");

  /**
   * Destructor drops the logger from the spd registry.
   */
  ~ConjunctionLogger() override;

  void writeConjunctions(size_t iteration,
                         const std::unordered_map<Particle *, std::tuple<Particle *, double>> &collisions) override;

 private:
  /**
   * Log the given arguments and the internal buffer to the csv file.
   * @param p1
   * @param p2
   */
  void log(size_t iteration, const Particle &p1, const Particle &p2, double distance);

  std::string _loggerName{"conjunctionLogger"};
};
