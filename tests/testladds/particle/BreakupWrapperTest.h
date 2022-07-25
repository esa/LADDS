/**
 * @file BreakupWrapperTest.h
 * @author F. Gratl
 * @date 24.01.22
 */

#pragma once
#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include "ladds/Simulation.h"
#include "ladds/TypeDefinitions.h"
#include "ladds/distributedMemParallelization/RegularGridDecomposition.h"
#include "ladds/io/Logger.h"

class BreakupWrapperTest : public testing::Test {
 public:
  BreakupWrapperTest() : logger("BreakupWrapperTestLogger"), simulation(logger) {
    logger.get()->set_level(LADDS::Logger::Level::off);

    // initialize a minimal default configuration
    config["autopas"]["cutoff"] = 3;
    config["sim"]["deltaT"] = 1.0;
    config["sim"]["maxAltitude"] = 85000.;
    config["sim"]["prop"]["coefficientOfDrag"] = 2.2;

    // optional parameters which are necessary for the tests here
    // no propagator components are activated (not even kepler) to keep trajectories simpler
    config["sim"]["iterations"] = 1;
    config["sim"]["breakup"]["enabled"] = true;

    // init decomp in cfg file correctly
    config["sim"]["decompositionType"] = "RegularGrid";

    configReader = std::make_unique<LADDS::ConfigReader>(config, logger);
    decomposition = std::make_unique<LADDS::RegularGridDecomposition>(*configReader);
    autopas = simulation.initAutoPas(*configReader, *decomposition);
  };

  YAML::Node config;

  LADDS::Logger logger;
  LADDS::Simulation simulation;
  std::unique_ptr<LADDS::ConfigReader> configReader;
  std::unique_ptr<LADDS::DomainDecomposition> decomposition;
  std::unique_ptr<AutoPas_t> autopas;
};
