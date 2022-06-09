/**
 * @file ConfigReader.cpp
 * @author F. Gratl
 * @date 20.12.21
 */

#include "ConfigReader.h"

#include <fstream>
#include <iomanip>

#include "ladds/io/hdf5/HDF5Definitions.h"

namespace LADDS {

ConfigReader::ConfigReader(const std::string &configPath, const Logger &logger, int rank, int numRanks)
    : config(loadConfig(configPath, logger)), logger(logger) {
  const auto constellationList = get<std::string>("io/constellationList", "", true);
  const auto numConstellations =
      constellationList.empty() ? 0 : std::count(constellationList.begin(), constellationList.end(), ';') + 1;
  const auto lengthIDRange = std::numeric_limits<HDF5Definitions::IntType>::max() / (numRanks + numConstellations);
  nextSafeParticleID = lengthIDRange * rank;
  lastSafeParticleID = lengthIDRange * (rank + 1) - 1;
}
size_t ConfigReader::newParticleID() {
  if (nextSafeParticleID > lastSafeParticleID) {
    throw std::runtime_error("No particle IDs left in this rank's ID space!");
  }
  return nextSafeParticleID++;
}

YAML::Node ConfigReader::loadConfig(const std::string &cfgFilePath, const Logger &logger) {
  YAML::Node file;
  try {
    file = YAML::LoadFile(cfgFilePath);
  } catch (YAML::Exception &e) {
    logger.log(Logger::Level::warn, "Error loading cfg, loading default file.\n{}", e.msg);
    try {
      file = YAML::LoadFile(defaultCfgPath);
    } catch (YAML::Exception &e) {
      logger.log(Logger::Level::err, "No default file file found. Should be cfg/default_cfg.yaml. Exiting...");
      exit(1);
    }
  }
  return file;
}

void ConfigReader::printParsedValues() {
  int maxKeyLength = 0;
  for (const auto &[key, value] : parsedValues) {
    maxKeyLength = std::max(static_cast<int>(key.size()), maxKeyLength);
  }

  std::stringstream ss{};
  for (const auto &[key, value] : parsedValues) {
    ss << "\n  " << std::left << std::setw(maxKeyLength + 2) << key << ": " << std::right << value;
  }
  SPDLOG_LOGGER_INFO(logger.get(), "Configuration:{}", ss.str());
}

bool ConfigReader::defines(const std::string &valuePath, bool suppressWarning) {
  std::vector<std::string> valuePathVec = autopas::utils::StringUtils::tokenize(valuePath, "/");
  YAML::Node node = Clone(config);
  for (const auto &dir : valuePathVec) {
    // as soon as anything in the path is undefined abort and either use the fallback value or abort.
    if (not node[dir].IsDefined()) {
      if (not suppressWarning) {
        SPDLOG_LOGGER_WARN(logger.get(), "Config value \"{}\" not defined!", valuePath);
      }
      return false;
    }
    node = node[dir];
  }
  return true;
}

void ConfigReader::dumpConfig(const std::string &filename) const {
  std::ofstream file(filename);
  file << config << std::endl;
  file.close();
}

const Logger &ConfigReader::getLogger() const {
  return logger;
}

size_t ConfigReader::getFirstIterationNr() {
  // either iteration of HDF5 checkpoint or 0
  return defines("io/hdf5", true) ? get<size_t>("io/hdf5/checkpoint/iteration", -1, true) + 1 : 0;
}

size_t ConfigReader::getLastIterationNr() {
  return get<size_t>("sim/iterations") + getFirstIterationNr() - 1;
}

}  // namespace LADDS