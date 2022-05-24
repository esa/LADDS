/**
 * @file RegularGridDecompositionLogger.h
 * @author F. Gratl
 * @date 24.05.22
 */

#pragma once

#include <autopas/selectors/Configuration.h>

#include <string>

#include "DecompositionLogger.h"
#include "ladds/io/ConfigReader.h"
#include "ladds/distributedMemParallelization/RegularGridDecomposition.h"

namespace LADDS {
class RegularGridDecompositionLogger : public DecompositionLogger {
 public:
  RegularGridDecompositionLogger(ConfigReader &config, const LADDS::RegularGridDecomposition &decomposition);

  void writeMetafile(size_t iteration) const override;

  void writePayload(size_t iteration, const autopas::Configuration &autoPasConfig) const override;

 private:
  std::string loggerName;

  /**
   * Reference to the domain decomposition.
   */
  const RegularGridDecomposition &decomposition;

  /**
   * Stores the maximum number of digits an iteration can have.
   * This is used to determine the number of leading zeros for each timestep record.
   */
  size_t maxDigitsIterations;

  /**
   * Generates a filename for the metadata file (pvts).
   * @param iteration
   * @return
   */
  std::string filenameMetadata(size_t iteration) const;

  /**
   * Generates a filename for the payload file (vts).
   * @param iteration
   * @param rank
   * @return
   */
  std::string filenamePayload(size_t iteration, int rank = -1) const;
};
}  // namespace LADDS