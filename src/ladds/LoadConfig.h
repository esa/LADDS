#pragma once

#include <yaml-cpp/yaml.h>

#include "Logger.h"

namespace LoadConfig {

constexpr auto defaultCfgPath = "../../cfg/default_cfg.yaml";

/**
 *  Loads the config file, falls back to default cfg if not found.
 */
[[nodiscard]] YAML::Node loadConfig(const std::string &cfgFilePath, const Logger &logger) {
  YAML::Node config;
  try {
    config = YAML::LoadFile(cfgFilePath);
  } catch (YAML::Exception &e) {
    std::cout << e.msg << std::endl;
    logger.log(Logger::Level::warn, "Error loading cfg, loading default config.");
    try {
      config = YAML::LoadFile(defaultCfgPath);
    } catch (YAML::Exception &e) {
      logger.log(Logger::Level::err, "No default config file found. Should be ../../cfg/default_cfg.yaml. Exiting...");
      exit(1);
    }
  }
  return config;
}
}  // namespace LoadConfig