/**
 * @file SatelliteLoader.cpp
 * @author F. Gratl
 * @date 06.12.21
 */

#include "SatelliteLoader.h"

#include "ConfigReader.h"
#include "DatasetReader.h"
#include "Logger.h"
#include "ladds/io/hdf5/HDF5Reader.h"

void SatelliteLoader::loadSatellites(AutoPas_t &autopas, ConfigReader &config, const Logger &logger) {
  std::vector<Particle> satellites;

  // load CSV ...
  if (const auto posFilePathCfg = config.get<std::string>("io/posFileName", ""),
      velFilePathCfg = config.get<std::string>("io/velFileName", "");
      not posFilePathCfg.empty() and not velFilePathCfg.empty()) {
    const auto posFilePath = std::string(DATADIR) + posFilePathCfg;
    const auto velFilePath = std::string(DATADIR) + velFilePathCfg;

    SPDLOG_LOGGER_INFO(
        logger.get(), "Loading scenario from CSV\nPositions: {}\nVelocities: {}", posFilePath, velFilePath);

    // Read in scenario
    satellites = DatasetReader::readDataset(posFilePath, velFilePath);
    SPDLOG_LOGGER_DEBUG(logger.get(), "Parsed {} satellites", satellites.size());
  } else
      // ... or load checkpoint ...
      if (const auto checkpointPathCfg = config.get<std::string>("io/checkpoint/file", "");
          not checkpointPathCfg.empty()) {
    const auto checkpointPath = std::string(DATADIR) + checkpointPathCfg;
    SPDLOG_LOGGER_INFO(logger.get(), "Loading scenario from HDF5 checkpoint\nFile: {}", checkpointPath);

    HDF5Reader hdfReader(checkpointPath);
    // either load the given iteration or fall back to the last iteration stored in the file
    auto iteration = config.get<size_t>("io/checkpoint/iteration", hdfReader.readLastIterationNr());
    satellites = hdfReader.readParticles(iteration);
  } else {
    // ... or fail
    throw std::runtime_error("No valid input option found! Exiting...");
  }

  // load particle vector into autopas while checking that they are within the desired altitude
  const auto maxAltitude = config.get<double>("sim/maxAltitude");
  const auto maxAltitudeSquared = maxAltitude * maxAltitude;
  double minAltitudeFound{std::numeric_limits<double>::max()};
  double maxAltitudeFound{0.};
  for (const auto &particle : satellites) {
    const auto &pos = particle.getPosition();
    const double altitudeSquared = autopas::utils::ArrayMath::dot(pos, pos);
    minAltitudeFound = std::min(minAltitudeFound, altitudeSquared);
    maxAltitudeFound = std::max(maxAltitudeFound, altitudeSquared);
    if (altitudeSquared < maxAltitudeSquared) {
      autopas.addParticle(particle);
    } else {
      SPDLOG_LOGGER_WARN(logger.get(),
                         "Particle NOT added because its altitudeSquared was too high!\n"
                         "Max allowed: {}\n"
                         "Actual: {}\n"
                         "{})",
                         maxAltitude,
                         sqrt(altitudeSquared),
                         particle.toString());
    }
  }
  minAltitudeFound = sqrt(minAltitudeFound);
  maxAltitudeFound = sqrt(maxAltitudeFound);
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