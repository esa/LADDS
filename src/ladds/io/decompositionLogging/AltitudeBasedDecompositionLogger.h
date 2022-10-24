/**
 * @file AltitudeBasedDecompositionLogger.h
 * @author F. Gratl
 * @date 24.10.22
 */

#pragma once

#include <string>

#include "DecompositionLoggerParametrized.h"
#include "ladds/distributedMemParallelization/AltitudeBasedDecomposition.h"
#include "ladds/io/ConfigReader.h"

namespace LADDS {
class AltitudeBasedDecompositionLogger : public DecompositionLoggerParametrized<AltitudeBasedDecomposition> {
 public:
  AltitudeBasedDecompositionLogger(ConfigReader &config, const LADDS::AltitudeBasedDecomposition &decomposition);

  void writeMetafile(size_t iteration) const override;

  void writePayload(size_t iteration, const autopas::Configuration &autoPasConfig) const override;
};
}  // namespace LADDS