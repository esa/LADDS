/**
 * @file SatelliteLoader.h
 * @author F. Gratl
 * @date 06.12.21
 */

#pragma once

#include "ladds/TypeDefinitions.h"
#include "ladds/particle/Constellation.h"
#include "yaml-cpp/yaml.h"
#include "Logger.h"

namespace SatelliteLoader {
/**
 * Load the particles from the input csv files in the config as particles into AutoPas.
 *
 * @note Paths for csv files are relative to ladds/data!
 *
 * @param autopas
 * @param config
 */
void loadSatellites(AutoPas_t &autopas, const YAML::Node &config, const Logger &logger);

/**
 * Parse constellation information and prepare satellites for insertion.
 * @param config
 * @return Vector of Constellations
 */
std::vector<Constellation> loadConstellations(const YAML::Node &config, const Logger &logger);

}  // namespace SatelliteLoader
