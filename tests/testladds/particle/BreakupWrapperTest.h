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
#include "ladds/io/Logger.h"

class BreakupWrapperTest : public testing::Test {
 public:
  BreakupWrapperTest() : logger("BreakupWrapperTestLogger"), simulation(logger) {
    logger.get()->set_level(Logger::Level::off);

    // initialize a minimal default configuration
    config["autopas"]["cutoff"] = 3;
    config["sim"]["deltaT"] = 1.0;
    config["sim"]["maxAltitude"] = 85000.;
    config["sim"]["coefficientOfDrag"] = 2.2;

    // optional parameters which are necessary for the tests here
    // no propagator components are activated (not even kepler) to keep trajectories simpler
    config["sim"]["iterations"] = 1;
    config["sim"]["conjunctionThreshold"] = 0.01;
    config["sim"]["breakup"]["enabled"] = true;

    configReader = std::make_unique<ConfigReader>(config, logger);
    autopas = simulation.initAutoPas(*configReader);
  };

  YAML::Node config;

  Logger logger;
  Simulation simulation;
  std::unique_ptr<ConfigReader> configReader;
  std::unique_ptr<AutoPas_t> autopas;
};
