/**
 * @file main.cpp
 * @author F. Gratl
 * @date 28.06.21
 */

#include <iostream>
#include <string>

#include "Simulation.h"
#include "ladds/Version.h"
#include "ladds/io/ConfigReader.h"
#include "ladds/io/Logger.h"
#include "autopas/utils/WrapMPI.h"

/**
 * The entrypoint of the program.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv) {
  LADDS::Logger logger;

  autopas::AutoPas_MPI_Init(&argc, &argv);

  // Default config path
  if (argc != 2) {
    logger.log(LADDS::Logger::Level::critical, "No config given!");
    return -1;
  }

  // Read in config
  const auto *cfgFilePath = argv[1];
  auto config = LADDS::ConfigReader(cfgFilePath, logger);

  logger.get()->set_level(spdlog::level::from_str(config.get<std::string>("sim/logLevel", "info")));
  SPDLOG_LOGGER_INFO(logger.get(), "LADDS version: {}", LADDS_VERSION);
  SPDLOG_LOGGER_INFO(logger.get(), "Config loaded.");

  LADDS::Simulation simulation(logger);

  simulation.run(config);

  autopas::AutoPas_MPI_Finalize();
  return 0;
}
