/**
 * @file ConfigReader.cpp
 * @author F. Gratl
 * @date 20.12.21
 */

#include "ConfigReader.h"

#include <iomanip>

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
