/**
 * @file DatasetReader.cpp
 * @author F. Gratl
 * @date 01.12.21
 */

#include "DatasetReader.h"

std::vector<Particle> DatasetReader::readDataset(const std::string &csvFilepath) {
  CSVReader<size_t,       // id
            std::string,  // name
            double,       // mass
            double,       // radius
            std::string,  // Particle::ActivityState
            double,       // r x
            double,       // r y
            double,       // r z
            double,       // v x
            double,       // v y
            double        // v z
            >
      csvReader{csvFilepath, true};
  const auto parsedData = csvReader.getLines();

  std::vector<Particle> particleCollection;
  particleCollection.reserve(parsedData.size());

  size_t particleId = 0;
  for (const auto &[id, name, mass, radius, activityState, rX, rY, rZ, vX, vY, vZ] : parsedData) {
    particleCollection.emplace_back(std::array<double, 3>{rX, rY, rZ}, std::array<double, 3>{vX, vY, vZ}, id);
  }

  return particleCollection;
}
