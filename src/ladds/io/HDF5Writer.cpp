/**
 * @file HDF5Writer.cpp
 * @author F. Gratl
 * @date 14.12.21
 */

#include "HDF5Writer.h"

void HDF5Writer::write(size_t iteration, const AutoPas_t &autopas) {
  std::vector<ParticleData> data;
  data.reserve(autopas.getNumberOfParticles());

  for (const auto &particle : autopas) {
    const auto &pos = particle.getR();
    const auto &vel = particle.getV();
    // pack data and make sure it is of the correct type
    data.emplace_back<ParticleData>({static_cast<decltype(ParticleData::rx)>(pos[0]),
                                     static_cast<decltype(ParticleData::ry)>(pos[1]),
                                     static_cast<decltype(ParticleData::rz)>(pos[2]),
                                     static_cast<decltype(ParticleData::vx)>(vel[0]),
                                     static_cast<decltype(ParticleData::vy)>(vel[1]),
                                     static_cast<decltype(ParticleData::vz)>(vel[2]),
                                     static_cast<decltype(ParticleData::id)>(particle.getID())});
  }

  const auto group = "Iteration" + std::to_string(iteration);
  _file.writeDataset(data, group + "/Particles", ParticleDataH5Type);
  _file.writeAttribute(iteration, "IterationNr", group);
}
