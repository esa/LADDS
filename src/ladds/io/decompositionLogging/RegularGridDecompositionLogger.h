/**
 * @file RegularGridDecompositionLogger.h
 * @author F. Gratl
 * @date 24.05.22
 */

#pragma once

#include <string>

#include "DecompositionLoggerParametrized.h"
#include "ladds/distributedMemParallelization/RegularGridDecomposition.h"
#include "ladds/io/ConfigReader.h"

namespace LADDS {
class RegularGridDecompositionLogger : public DecompositionLoggerParametrized<RegularGridDecomposition> {
 public:
  RegularGridDecompositionLogger(ConfigReader &config, const LADDS::RegularGridDecomposition &decomposition);

  void writeMetafile(size_t iteration) const override;

  void writePayload(size_t iteration, const autopas::Configuration &autoPasConfig) const override;
};
}  // namespace LADDS