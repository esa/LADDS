/**
 * @file SatelliteLoader.cpp
 * @author F. Gratl
 * @date 06.12.21
 */

#include "SatelliteLoader.h"

#include "DatasetReader.h"
#include "Logger.h"

void SatelliteLoader::loadSatellites(AutoPas_t &autopas, const YAML::Node &config, const Logger &logger) {
  // Read in scenario
  auto actualSatellites =
      DatasetReader::readDataset(std::string(DATADIR) + config["io"]["posFileName"].as<std::string>(),
                                 std::string(DATADIR) + config["io"]["velFileName"].as<std::string>());
  SPDLOG_LOGGER_DEBUG(logger.get(), "Parsed {} satellites", actualSatellites.size());

  const auto maxAltitude = config["sim"]["maxAltitude"].as<double>();
  double minAltitudeFound{std::numeric_limits<double>::max()};
  double maxAltitudeFound{0.};
  // Convert satellites to particles
  for (const auto &particle : actualSatellites) {
    auto pos = particle.getPosition();
    double altitude = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
    minAltitudeFound = std::min(minAltitudeFound, altitude);
    maxAltitudeFound = std::max(maxAltitudeFound, altitude);
    if (altitude < maxAltitude) {
      autopas.addParticle(particle);
    }
  }
  SPDLOG_LOGGER_INFO(logger.get(), "Min altitude is {}", minAltitudeFound);
  SPDLOG_LOGGER_INFO(logger.get(), "Max altitude is {}", maxAltitudeFound);
  SPDLOG_LOGGER_INFO(logger.get(), "Number of particles: {}", autopas.getNumberOfParticles());
}

std::vector<Constellation> SatelliteLoader::loadConstellations(const YAML::Node &config, const Logger &logger) {
  std::vector<Constellation> constellations;

  if (config["io"]["constellationList"].IsDefined()) {
    const auto insertionFrequency =
        config["io"]["constellationFrequency"].IsDefined() ? config["io"]["constellationFrequency"].as<int>() : 1;
    auto constellationDataStr = config["io"]["constellationList"].as<std::string>();
    // count constellation by counting ';'
    int nConstellations = 1;
    for (char con : constellationDataStr) {
      if (con == ';') {
        nConstellations++;
      }
    }

    // parse constellation info
    constellations.reserve(nConstellations);
    for (int i = 0; i < nConstellations; ++i) {
      unsigned long offset = constellationDataStr.find(';', 0);
      if (offset == 0) {
        constellations.emplace_back(Constellation(constellationDataStr, insertionFrequency));
        break;
      } else {
        constellations.emplace_back(Constellation(constellationDataStr.substr(0, offset), insertionFrequency));
        constellationDataStr.erase(0, offset + 1);
      }
    }

    size_t constellationTotalNumSatellites = 0;
    for (const auto &constellation : constellations) {
      constellationTotalNumSatellites += constellation.getConstellationSize();
    }

    SPDLOG_LOGGER_INFO(logger.get(),
                       "{} more particles will be added from {} constellations",
                       constellationTotalNumSatellites,
                       nConstellations);
  }
  return constellations;
}