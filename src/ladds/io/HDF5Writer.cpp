/**
 * @file HDF5Writer.cpp
 * @author F. Gratl
 * @date 14.12.21
 */

#include "HDF5Writer.h"

std::array<float, 3> toFloat(const std::array<double, 3> &arrDouble) {
  std::array<float, 3> arrFloat{};
  std::copy(arrDouble.begin(), arrDouble.end(), arrFloat.begin());
  return arrFloat;
}

void HDF5Writer::write(size_t iteration, const AutoPas_t &autopas) {
  //  std::vector<data_t> data;
  //  data.reserve(autopas.getNumberOfParticles());
  std::vector<std::array<float, 3>> positions;
  std::vector<std::array<float, 3>> velocities;
  std::vector<size_t> ids;
  positions.reserve(autopas.getNumberOfParticles());
  velocities.reserve(autopas.getNumberOfParticles());
  ids.reserve(autopas.getNumberOfParticles());

  for (const auto &particle : autopas) {
    //    data.emplace_back<data_t>({particle.getID(), particle.getR(), particle.getV()});
    positions.push_back(toFloat(particle.getR()));
    velocities.push_back(toFloat(particle.getV()));
    ids.push_back(particle.getID());
  }

  auto group = _file.createGroup("Iteration" + std::to_string(iteration));

  auto datasetPos = group.createDataSet<float>("Positions", HighFive::DataSpace(autopas.getNumberOfParticles(), 3));
  datasetPos.write(positions.data());
  auto datasetV = group.createDataSet<float>("Velocities", HighFive::DataSpace(autopas.getNumberOfParticles(), 3));
  datasetV.write(velocities.data());
  auto datasetId = group.createDataSet("IDs", ids);

  group.createAttribute("IterationNr", iteration);
}
