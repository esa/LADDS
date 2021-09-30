#pragma once

#include <breakupModel/input/CSVReader.h>

#include "Particle.h"

namespace DatasetReader {

/**
 * Reads the passed position and velocity csv files. Returns a vector of particles.
 */
std::vector<Particle> readDataset(const std::string position_filepath, const std::string velocity_filepath) {
  CSVReader<double, double, double> pos_csvReader{position_filepath, false};
  CSVReader<double, double, double> vel_csvReader{velocity_filepath, false};
  std::vector<Particle> particleCollection;

  auto positions = pos_csvReader.getLines();
  auto velocities = vel_csvReader.getLines();

  if (positions.size() != velocities.size()) {
    std::cout << "Error: Position and velocity file have different number of lines." << std::endl;
    return particleCollection;
  }

  for (unsigned i = 0; i < positions.size(); i++) {
    const auto &[x, y, z] = positions[i];
    const auto &[vx, vy, vz] = velocities[i];

    const std::array<double, 3> pos = {x, y, z};
    const std::array<double, 3> vel = {vx, vy, vz};
    particleCollection.push_back(Particle(pos, vel, i));
  }

  return particleCollection;
}
}  // namespace DatasetReader