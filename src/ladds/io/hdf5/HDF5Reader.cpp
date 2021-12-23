/**
 * @file HDF5Reader.cpp
 * @author F. Gratl
 * @date 23.12.21
 */

#include "HDF5Reader.h"

#include "HDF5Writer.h"
HDF5Reader::HDF5Reader(const std::string &filename)
#ifdef LADDS_HDF5
    : file(filename, h5pp::FilePermission::READONLY)
#endif
{
}

std::vector<Particle> HDF5Reader::readParticles(size_t iteration) {
  std::vector<Particle> particles{};
#ifdef LADDS_HDF5
  const auto pos = file.readDataset<std::vector<HDF5Writer::Vec3<HDF5Writer::FloatType>>>(
      HDF5Writer::groupParticleData + std::to_string(iteration) + HDF5Writer::datasetParticlePositions);
  const auto vel = file.readDataset<std::vector<HDF5Writer::Vec3<HDF5Writer::FloatType>>>(
      HDF5Writer::groupParticleData + std::to_string(iteration) + HDF5Writer::datasetParticleVelocities);
  const auto ids = file.readDataset<std::vector<HDF5Writer::IntType>>(
      HDF5Writer::groupParticleData + std::to_string(iteration) + HDF5Writer::datasetParticleIDs);

  particles.reserve(pos.size());
  for (size_t i = 0; i < pos.size(); ++i) {
    particles.emplace_back(std::array<double, 3>{pos[i].x, pos[i].y, pos[i].z},
                           std::array<double, 3>{vel[i].x, vel[i].y, vel[i].z},
                           ids[i]);
  }
#endif
  return particles;
}

std::vector<HDF5Writer::CollisionInfo> HDF5Reader::readCollisions(size_t iteration) {
  std::vector<HDF5Writer::CollisionInfo> collisions{};
#ifdef LADDS_HDF5
  file.readDataset(collisions,
                   HDF5Writer::groupCollisionData + std::to_string(iteration) + HDF5Writer::datasetCollisions);
#endif
  return collisions;
}
