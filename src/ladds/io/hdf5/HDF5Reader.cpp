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
#ifndef LADDS_HDF5
  throw std::runtime_error("LADDS was compiled without HDF5 support, so the HDF5Reader can't do anything!");
#endif
}

std::vector<Particle> HDF5Reader::readParticles(size_t iteration) const {
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

std::vector<HDF5Writer::CollisionInfo> HDF5Reader::readCollisions(size_t iteration) const {
  std::vector<HDF5Writer::CollisionInfo> collisions{};
#ifdef LADDS_HDF5
  file.readDataset(collisions,
                   HDF5Writer::groupCollisionData + std::to_string(iteration) + HDF5Writer::datasetCollisions);
#endif
  return collisions;
}

size_t HDF5Reader::readLastIterationNr() {
  std::vector<size_t> allIterations;
#ifdef LADDS_HDF5
  const auto groupStrLength = std::string(HDF5Writer::groupParticleData).size();
  // gather all ID datasets. They contain the iteration number in their path
  auto allIdDatasets = file.findDatasets(HDF5Writer::datasetParticleIDs);
  allIterations.reserve(allIdDatasets.size());
  // extract all iteration numbers as size_t
  std::transform(allIdDatasets.begin(), allIdDatasets.end(), std::back_inserter(allIterations), [&](auto &datasetName) {
    const auto posSecondSlash = datasetName.find('/', groupStrLength);
    const auto iterationStr = datasetName.substr(groupStrLength, posSecondSlash - groupStrLength);
    return strtoul(iterationStr.c_str(), nullptr, 10);
  });
  // sort numerically and return the highest
  std::sort(allIterations.begin(), allIterations.end());
#endif
  return allIterations.back();
}
