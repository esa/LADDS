/**
 * @file DecompositionLogger.h
 * @author F. Gratl
 * @date 24.05.22
 */

#pragma once

#include <utility>

#include "ladds/io/ConfigReader.h"

namespace LADDS {

class DecompositionLogger {
 public:
  DecompositionLogger(std::string loggerName, LADDS::ConfigReader &config, std::string fileExtension);

  virtual ~DecompositionLogger() = default;
  /**
   * Creates the pvts file that references all vts files
   * @note Should only be called by one rank!
   * @param iteration
   */
  virtual void writeMetafile(size_t iteration) const = 0;

  /**
   * Creates all vts files which contain the actual decomposition.
   * @param iteration
   * @param autoPasConfig
   */
  virtual void writePayload(size_t iteration, const autopas::Configuration &autoPasConfig) const = 0;

 protected:
  /**
   * Name for the underlying SPD logger.
   */
  std::string loggerName;

  /**
   * Stores the maximum number of digits an iteration can have.
   * This is used to determine the number of leading zeros for each timestep record.
   */
  size_t maxDigitsIterations;

  /**
   * Extension for output file names (e.g. vtk / vtu / vtp...).
   */
  std::string fileExtension;
};
}  // namespace LADDS