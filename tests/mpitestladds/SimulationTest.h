/**
 * @file SimulationTest.h
 * @author F. Gratl
 * @date 02.06.22
 */

#pragma once

#include <gtest/gtest.h>
#include <yaml-cpp/node/node.h>

#include "autopas/utils/WrapOpenMP.h"
#include "ladds/Simulation.h"
#include "ladds/distributedMemParallelization/DomainDecomposition.h"
#include "ladds/io/ConfigReader.h"
#include "ladds/io/Logger.h"

class SimulationTest : public testing::Test {
 public:
  SimulationTest()
      : maxThreadsBefore(autopas::autopas_get_max_threads()), logger("SimulationTestLogger"), simulation(logger) {
    // make sure to only use one thread
    autopas::autopas_set_num_threads(1);

    logger.get()->set_level(LADDS::Logger::Level::err);

    // initialize a minimal default configuration
    config["autopas"]["cutoff"] = 80.;
    config["sim"]["breakup"]["enabled"] = false;
    config["sim"]["deltaT"] = 1.0;
    config["sim"]["maxAltitude"] = 85000.;
    config["sim"]["prop"]["coefficientOfDrag"] = 2.2;

    // optional parameters which are necessary for the tests here
    config["io"]["constellationCutoff"] = constellationCutoff;
    config["sim"]["collisionDistanceFactor"] = 1.;
    config["sim"]["iterations"] = 1;
    config["sim"]["minAltitude"] = 150.;
    config["sim"]["prop"]["useKEPComponent"] = true;

    configReader = std::make_unique<LADDS::ConfigReader>(config, logger);

    decomposition = std::make_unique<LADDS::RegularGridDecomposition>(*configReader);

    autopas = simulation.initAutoPas(*configReader, *decomposition);
  }

  virtual ~SimulationTest() {
    // reset omp max threads
    autopas::autopas_set_num_threads(maxThreadsBefore);
  }

 public:
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
