/**
 * @file main.cpp
 * @author F. Gratl
 * @date 28.06.21
 */

#include <iostream>

#include "Simulation.h"
#include "ladds/io/LoadConfig.h"
#include "ladds/io/Logger.h"

/**
 * The entrypoint of the program.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv) {
  Logger logger;

  // Default config path
  if (argc != 2) {
    logger.log(Logger::Level::critical, "No config given!");
    return -1;
  }

  // Read in config
  const auto *cfgFilePath = argv[1];
  const auto config = LoadConfig::loadConfig(cfgFilePath, logger);

  logger.get()->set_level(spdlog::level::from_str(config["sim"]["logLevel"].as<std::string>()));
  SPDLOG_LOGGER_INFO(logger.get(), "Config loaded.");

  Simulation simulation(logger);

  simulation.run(config);

  return 0;
}
