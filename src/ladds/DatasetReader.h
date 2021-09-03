#pragma once

#include <breakupModel/input/CSVReader.h>

#include "Particle.h"

namespace DatasetReader {
std::vector<Particle> readDataset(const std::string position_filepath, const std::string velocity_filepath) {
  std::cout << position_filepath << std::endl;

  CSVReader pos_csvReader{position_filepath, false};
  CSVReader vel_csvReader{velocity_filepath, false};
  std::vector<Particle> particleCollection;

  auto positions = pos_csvReader.getLines();
  auto velocities = vel_csvReader.getLines();

  for (const auto line : positions) {
    std::cout << std::get<0>(line) << std::endl;
  }

  return particleCollection;
}
}  // namespace DatasetReader