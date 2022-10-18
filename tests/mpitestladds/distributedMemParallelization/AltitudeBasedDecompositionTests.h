/**
 * @file AltitudeBasedDecompositionTests.h
 * @author P. Gomez
 * @date 2022-08-01
 */

#pragma once

#include <gtest/gtest.h>

#include "ladds/Simulation.h"
#include "ladds/distributedMemParallelization/AltitudeBasedDecomposition.h"
#include "ladds/io/ConfigReader.h"
#include "ladds/io/Logger.h"

class AltitudeBasedDecompositionTests : public testing::Test {
 public:
  AltitudeBasedDecompositionTests();

  ~AltitudeBasedDecompositionTests() override;

  int maxThreadsBefore;
  YAML::Node config;

  LADDS::Logger logger;
  LADDS::Simulation simulation;
  std::unique_ptr<LADDS::ConfigReader> configReader;
  std::unique_ptr<AutoPas_t> autopas;
  std::unique_ptr<LADDS::AltitudeBasedDecomposition> decomposition;
  static constexpr std::array<double, 3> zeroVec{0., 0., 0.};
  static constexpr double constellationCutoff = 0.04;
};
