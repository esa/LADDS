/**
 * @file SatelliteLoader.cpp
 * @author F. Gratl
 * @date 06.12.21
 */

#include "SatelliteLoader.h"

#include "ConfigReader.h"
#include "DatasetReader.h"
#include "Logger.h"

void SatelliteLoader::loadSatellites(AutoPas_t &autopas, ConfigReader &config, const Logger &logger) {
  // Read in scenario
  auto actualSatellites = DatasetReader::readDataset(std::string(DATADIR) + config.get<std::string>("io/posFileName"),
                                                     std::string(DATADIR) + config.get<std::string>("io/velFileName"));
  SPDLOG_LOGGER_DEBUG(logger.get(), "Parsed {} satellites", actualSatellites.size());

  const auto maxAltitude = config.get<double>("sim/maxAltitude");
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

std::vector<Constellation> SatelliteLoader::loadConstellations(ConfigReader &config, const Logger &logger) {
  std::vector<Constellation> constellations;
  auto constellationList = config.get<std::string>("io/constellationList", "", true);
  // altitudeSpread = 3 * sigma , altitudeDeviation = sigma (= standardDeviation)
  auto altitudeDeviation = config.get<double>("io/altitudeSpread", 0.0) / 3.0;
  if (not constellationList.empty()) {
    const auto insertionFrequency = config.get<int>("io/constellationFrequency", 1);
    auto constellationDataStr = config.get<std::string>("io/constellationList");
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
        constellations.emplace_back(Constellation(constellationDataStr, insertionFrequency, altitudeDeviation));
        break;
      } else {
        constellations.emplace_back(
            Constellation(constellationDataStr.substr(0, offset), insertionFrequency, altitudeDeviation));
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