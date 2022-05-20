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

namespace LADDS {

void SatelliteLoader::loadSatellites(AutoPas_t &autopas, ConfigReader &config, const Logger &logger) {
  std::vector<Particle> satellites;

  // if the value is not set this is the only place that should define the default value as it is the first to access it
  const auto coefficientOfDrag = config.get<double>("sim/prop/coefficientOfDrag", 2.2);

  const auto csvFileName = config.get<std::string>("io/csv/fileName", "", true);
  // avoid creating io/hdf5 in the config if it doesn't exist
  const auto checkpointPathCfg =
      config.defines("io/hdf5", true) ? config.get<std::string>("io/hdf5/checkpoint/file", "", true) : "";
  // sanity checks
  if (csvFileName.empty() and checkpointPathCfg.empty()) {
    throw std::runtime_error("No valid input option found! Exiting...");
  }
  if (not csvFileName.empty() and not checkpointPathCfg.empty()) {
    throw std::runtime_error("Ambiguous input option found! Both CSV and HDF5 checkpoint found. Exiting...");
  }

  // load CSV ...
  if (not csvFileName.empty()) {
    // TODO: Remove DATADIR functionality
    const auto csvFilePath = std::string(DATADIR) + csvFileName;

    SPDLOG_LOGGER_INFO(logger.get(), "Loading scenario from CSV: {}", csvFilePath);

    // Read in scenario
    satellites = DatasetReader::readDataset(csvFilePath, coefficientOfDrag);
    SPDLOG_LOGGER_DEBUG(logger.get(), "Parsed {} satellites", satellites.size());
  }
  // ... or load checkpoint ...
  if (not checkpointPathCfg.empty()) {
    // TODO: Remove DATADIR functionality
    const auto checkpointPath = std::string(DATADIR) + checkpointPathCfg;
    SPDLOG_LOGGER_INFO(logger.get(), "Loading scenario from HDF5 checkpoint\nFile: {}", checkpointPath);

    HDF5Reader hdfReader(checkpointPath);
    // either load the given iteration or fall back to the last iteration stored in the file
    auto iteration = config.get<size_t>("io/hdf5/checkpoint/iteration", hdfReader.readLastIterationNr());
    satellites = hdfReader.readParticles(iteration, coefficientOfDrag);
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
  if (not constellationList.empty()) {
    const auto insertionFrequency = config.get<int>("io/constellationFrequency", 1);
    auto constellationDataStr = config.get<std::string>("io/constellationList");
    // count constellation by counting ';'
    int nConstellations = std::count(constellationDataStr.begin(), constellationDataStr.end(), ';') + 1;

    // parse constellation info
    constellations.reserve(nConstellations);
    for (int i = 0; i < nConstellations; ++i) {
      unsigned long offset =
          (i == nConstellations - 1) ? constellationDataStr.size() : constellationDataStr.find(';', 0);
      std::string constellationDir = constellationDataStr.substr(0, offset);

      ConfigReader constellationConfig =
          ConfigReader(std::string(DATADIR) + constellationDir + "/shells_" + constellationDir + ".yaml", logger);

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
                       "{} more particles will be added from {} constellations",
                       constellationTotalNumSatellites,
                       nConstellations);

    for (const auto &c : constellations) {
      SPDLOG_LOGGER_INFO(
          logger.get(),
          "{}: insertion starts at iteration: {}, is fully deployed within {} iterations, inserts {} satellites",
          c.getConstellationName(),
          c.getStartTime(),
          c.getDuration(),
          c.getConstellationSize());
    }
  }
  return constellations;
}

}  // namespace LADDS