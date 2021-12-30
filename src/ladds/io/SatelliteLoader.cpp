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
  auto constellationList =
      config["io"]["constellationList"].IsDefined() ? config["io"]["constellationList"].as<std::string>() : "";
  if (!constellationList.empty()) {
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
      unsigned long offset =
          (i == nConstellations - 1) ? constellationDataStr.size() : constellationDataStr.find(';', 0);
      std::string constellationDir = constellationDataStr.substr(0, offset);

      YAML::Node constellationConfig;
      try {
        constellationConfig =
            YAML::LoadFile(std::string(DATADIR) + constellationDir + "/shells_" + constellationDir + ".yaml");
      } catch (YAML::Exception &e) {
        std::cout << e.msg << std::endl;
        logger.log(Logger::Level::warn, "Error loading cfg, Exiting...");
        exit(1);
      }

      constellations.emplace_back(Constellation(constellationConfig, config));
      if (i != nConstellations - 1) {
        constellationDataStr.erase(0, offset + 1);
      }
    }

    size_t constellationTotalNumSatellites = 0;
    for (const auto &constellation : constellations) {
      constellationTotalNumSatellites += constellation.getConstellationSize();
    }

    SPDLOG_LOGGER_INFO(logger.get(),
                       "{} more satellites will be added from {} constellations:",
                       constellationTotalNumSatellites,
                       nConstellations);
  }
  for(auto & c : constellations){
      SPDLOG_LOGGER_INFO(logger.get(),
                         "{}: insertion starts at iteration: {}, is fully deployed within {} iterations, inserts {} satellites",
                         c.getConstellationName(),
                         c.getStartTime(),
                         c.getDuration(),
                         c.getConstellationSize());
  }
  return constellations;
}