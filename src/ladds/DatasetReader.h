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

  particleCollection.reserve(positions.size())

      size_t particleId = 0;
  std::transform(positions.begin(), positions.end(), velocities.begin(), std::back_insert_iterator(particleCollection),
                 [&](const auto &pos, const auto &vel) { return Particle(pos, vel, particleId++); });
  return particleCollection;
}
}  // namespace DatasetReader