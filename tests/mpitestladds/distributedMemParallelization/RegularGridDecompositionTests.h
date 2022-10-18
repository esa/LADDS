/**
 * @file RegularGridDecompositionTests.h
 * @author F. Gratl
 * @date 02.06.22
 */

#pragma once

#include <gtest/gtest.h>
#include <yaml-cpp/node/node.h>

#include "ladds/Simulation.h"
#include "ladds/distributedMemParallelization/DomainDecomposition.h"
#include "ladds/io/ConfigReader.h"
#include "ladds/io/Logger.h"

class RegularGridDecompositionTests : public testing::Test {
 public:
  RegularGridDecompositionTests();

  ~RegularGridDecompositionTests() override;

  int maxThreadsBefore;
  YAML::Node config;

  LADDS::Logger logger;
  LADDS::Simulation simulation;
  std::unique_ptr<LADDS::ConfigReader> configReader;
  std::unique_ptr<AutoPas_t> autopas;
  std::unique_ptr<LADDS::RegularGridDecomposition> decomposition;
  static constexpr std::array<double, 3> zeroVec{0., 0., 0.};
  static constexpr double constellationCutoff = 0.04;
};
