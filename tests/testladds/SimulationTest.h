/**
 * @file SimulationTest.h
 * @author F. Gratl
 * @date 03.12.21
 */

#pragma once

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include "ladds/Simulation.h"
#include "ladds/TypeDefinitions.h"
#include "ladds/io/Logger.h"

/**
 * position, positionIsSafe
 */
using ParameterTuple = std::tuple<std::array<double, 3>, bool>;

class SimulationTest : public testing::TestWithParam<ParameterTuple> {
 public:
  SimulationTest() : logger("SimulationTestLogger"), simulation(logger) {
    logger.get()->set_level(Logger::Level::off);

    // initialize a minimal default configuration
    config["sim"]["maxAltitude"] = 85000.;
    config["sim"]["deltaT"] = 1.0;
    config["autopas"]["cutoff"] = 0.02;
    config["autopas"]["skin"] = 0.2;
    config["autopas"]["rebuildFrequency"] = 20;
    config["autopas"]["desiredCellsPerDimension"] = 32;
    config["io"]["output_file"] = "test";
    config["io"]["constellationCutoff"] = constellationCutoff;
    config["sim"]["prop"]["useKEPComponent"] = true;

    configReader = std::make_unique<ConfigReader>(config, logger);
    autopas = simulation.initAutoPas(*configReader);
  }

  /**
   * Helper struct for pretty printing the generated test names.
   */
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      const auto &[pos, posIsSafe] = static_cast<ParamType>(info.param);

      std::stringstream ss;
      ss << "pos_" << autopas::utils::ArrayUtils::to_string(pos, "x", {"", ""}) << "_safe_" << std::boolalpha
         << posIsSafe;
      auto str = ss.str();
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      std::replace(str.begin(), str.end(), ' ', '_');
      return str;
    }
  };

  YAML::Node config;

  Logger logger;
  Simulation simulation;
  std::unique_ptr<ConfigReader> configReader;
  std::unique_ptr<AutoPas_t> autopas;
  static constexpr std::array<double, 3> zeroVec{0., 0., 0.};
  static constexpr std::array<double, 3> testCheckedInsertParticlePos{6871, 0, 0};
  static constexpr double constellationCutoff = 0.04;
};
