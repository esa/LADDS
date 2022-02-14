/**
 * @file DatasetReader.cpp
 * @author F. Gratl
 * @date 01.12.21
 */

#include "DatasetReader.h"

#include <satellitePropagator/physics/AccelerationComponents/DragComponent.h>

std::vector<Particle> DatasetReader::readDataset(const std::string &csvFilepath, double coefficientOfDrag) {
  CSVReader<size_t,                   // id
            std::string,              // cosparId
            std::string,              // name
            double,                   // bstar
            double,                   // mass
            double,                   // radius
            Particle::ActivityState,  // Particle::ActivityState
            double,                   // r x
            double,                   // r y
            double,                   // r z
            double,                   // v x
            double,                   // v y
            double                    // v z
            >
      csvReader{csvFilepath, true};
  const auto parsedData = csvReader.getLines();

  std::vector<Particle> particleCollection;
  particleCollection.reserve(parsedData.size());

  for (const auto &[id, cosparId, name, bstar, mass, radius, activityState, rX, rY, rZ, vX, vY, vZ] : parsedData) {
    const double bcInv = Particle::calculateBcInv(bstar, mass, radius, coefficientOfDrag);
    particleCollection.emplace_back(std::array<double, 3>{rX, rY, rZ},
                                    std::array<double, 3>{vX, vY, vZ},
                                    id,
                                    cosparId,
                                    activityState,
                                    mass,
                                    radius,
                                    bcInv);
  }

  return particleCollection;
}
