/**
 * @file DatasetReader.cpp
 * @author F. Gratl
 * @date 01.12.21
 */

#include "DatasetReader.h"

std::vector<Particle> DatasetReader::readDataset(const std::string &position_filepath,
                                                 const std::string &velocity_filepath) {
  CSVReader<double, double, double> pos_csvReader{position_filepath, false};
  CSVReader<double, double, double> vel_csvReader{velocity_filepath, false};
  std::vector<Particle> particleCollection;

  auto positions = pos_csvReader.getLines();
  auto velocities = vel_csvReader.getLines();

  if (positions.size() != velocities.size()) {
    std::cout << "Error: Position and velocity file have different number of lines." << std::endl;
    return particleCollection;
  }

  particleCollection.reserve(positions.size());

  size_t particleId = 0;
  std::transform(positions.begin(),
                 positions.end(),
                 velocities.begin(),
                 std::back_insert_iterator<std::vector<Particle>>(particleCollection),
                 [&](const auto &pos, const auto &vel) {
                   const auto &[x, y, z] = pos;
                   const auto &[vx, vy, vz] = vel;

                   const std::array<double, 3> posArray = {x, y, z};
                   const std::array<double, 3> velArray = {vx, vy, vz};
                   // TODO: parse activity information and set it here!
                   return Particle(posArray, velArray, particleId++, Particle::ActivityState::passive);
                 });
  return particleCollection;
}
