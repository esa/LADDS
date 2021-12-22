/**
 * @file SimulationTest.h
 * @author F. Gratl
 * @date 03.12.21
 */

#pragma once

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

class SimulationTest : public testing::Test {
public:
    SimulationTest() {
        // initialize a minimal default configuration
        config["sim"]["maxAltitude"] = 85000.;
        config["sim"]["deltaT"] = 1.0;
        config["autopas"]["cutoff"] = 0.02;
        config["autopas"]["skin"] = 0.2;
        config["autopas"]["rebuildFrequency"] = 20;
        config["autopas"]["desiredCellsPerDimension"] = 32;
        config["io"]["output_file"] = "test";
        config["io"]["constellationCutoff"] = 0.04;
        config["io"]["posFileName"] = "pos_test.csv";
        config["io"]["velFileName"] = "v_test.csv";
        config["sim"]["prop"]["useKEPComponent"] = true;
    }

    YAML::Node config;
};
