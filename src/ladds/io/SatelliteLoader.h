/**
 * @file SatelliteLoader.h
 * @author F. Gratl
 * @date 06.12.21
 */

#pragma once

#include "ConfigReader.h"
#include "Logger.h"
#include "ladds/TypeDefinitions.h"
#include "ladds/distributedMemParallelization/DomainDecomposition.h"
#include "ladds/particle/Constellation.h"
#include "yaml-cpp/yaml.h"

namespace LADDS::SatelliteLoader {
/**
 * Load the particles from the input csv files in the config as particles into AutoPas.
 *
 * @note Paths for csv files are relative to ladds/data!
 *
 * @param autopas
 * @param config
 * @param decomp Needed for the communicator
 */
void loadSatellites(AutoPas_t &autopas, ConfigReader &config, DomainDecomposition &decomp);

/**
 * Adds the passed particles to the autopas container.
 *
 * @param autopas
 * @param satellites
 * @param decomp
 * @param config
 * */
void addSatellitesToAutoPas(AutoPas_t &autopas,
                            std::vector<Particle> &satellites,
                            DomainDecomposition &decomp,
                            ConfigReader &config);

/**
 * Parse constellation information and prepare satellites for insertion.
 * @param config
 * @return Vector of Constellations
 */
std::vector<Constellation> loadConstellations(ConfigReader &config, const Logger &logger);

}  // namespace LADDS::SatelliteLoader
