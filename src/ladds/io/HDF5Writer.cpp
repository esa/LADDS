/**
 * @file HDF5Writer.cpp
 * @author F. Gratl
 * @date 14.12.21
 */

#include "HDF5Writer.h"

void HDF5Writer::write(size_t iteration, const AutoPas_t &autopas) {
  struct point3D {
    double x, y, z;
  };

  //  std::vector<data_t> data;
  //  data.reserve(autopas.getNumberOfParticles());
  std::vector<point3D> positions;
  std::vector<point3D> velocities;
  std::vector<size_t> ids;
  positions.reserve(autopas.getNumberOfParticles());
  velocities.reserve(autopas.getNumberOfParticles());
  ids.reserve(autopas.getNumberOfParticles());

  for (const auto &particle : autopas) {
    const auto &pos = particle.getR();
    positions.emplace_back<point3D>({pos[0], pos[1], pos[2]});
    const auto &vel = particle.getV();
    velocities.emplace_back<point3D>({vel[0], vel[1], vel[2]});
    ids.push_back(particle.getID());
  }

  const auto group = "Iteration" + std::to_string(iteration);
  _file.writeDataset(positions, group + "/Positions");
  _file.writeDataset(velocities, group + "/Velocities");
  _file.writeDataset(ids, group + "/IDs");
  _file.writeAttribute(iteration, "IterationNr", group);
}
