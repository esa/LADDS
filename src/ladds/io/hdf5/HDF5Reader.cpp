/**
 * @file HDF5Reader.cpp
 * @author F. Gratl
 * @date 23.12.21
 */

#include "HDF5Reader.h"

#include "HDF5Definitions.h"

HDF5Reader::HDF5Reader(const std::string &filename)
#ifdef LADDS_HDF5
    : file(filename, h5pp::FilePermission::READONLY)
#endif
{
#ifndef LADDS_HDF5
  throw std::runtime_error("LADDS was compiled without HDF5 support, so the HDF5Reader can't do anything!");
#endif
}

std::vector<Particle> HDF5Reader::readParticles(size_t iteration, double coefficientOfDrag) const {
  std::vector<Particle> particles{};
#ifdef LADDS_HDF5
  const auto pos = file.readDataset<std::vector<HDF5Definitions::Vec3<HDF5Definitions::FloatType>>>(
      HDF5Definitions::groupParticleData + std::to_string(iteration) + HDF5Definitions::datasetParticlePositions);
  const auto vel = file.readDataset<std::vector<HDF5Definitions::Vec3<HDF5Definitions::FloatType>>>(
      HDF5Definitions::groupParticleData + std::to_string(iteration) + HDF5Definitions::datasetParticleVelocities);
  const auto ids = file.readDataset<std::vector<HDF5Definitions::IntType>>(
      HDF5Definitions::groupParticleData + std::to_string(iteration) + HDF5Definitions::datasetParticleIDs);
  const auto pathParticleConstantProperties =
      std::string(HDF5Definitions::groupParticleData) + HDF5Definitions::tableParticleConstantProperties;
  const auto constantPropertiesVec =
      file.readTableRecords<std::vector<HDF5Definitions::ParticleConstantProperties>>(pathParticleConstantProperties);

  // convert static data to a map for easier access
  const auto constantPropertiesMap = [&]() {
    // id is redundant in this datastructure but easier maintainable this way
    std::unordered_map<size_t, HDF5Definitions::ParticleConstantProperties> map;
    map.reserve(constantPropertiesVec.size());
    for (const auto &data : constantPropertiesVec) {
      map[data.id] = data;
    }
    return map;
  }();

  particles.reserve(pos.size());
  for (size_t i = 0; i < pos.size(); ++i) {
    const auto &constantProperties = constantPropertiesMap.at(ids[i]);
    particles.emplace_back(std::array<double, 3>{pos[i].x, pos[i].y, pos[i].z},
                           std::array<double, 3>{vel[i].x, vel[i].y, vel[i].z},
                           constantProperties.id,
                           static_cast<Particle::ActivityState>(constantProperties.activityState),
                           constantProperties.mass,
                           constantProperties.radius,
                           coefficientOfDrag);
  }
#endif
  return particles;
}

std::vector<HDF5Definitions::CollisionInfo> HDF5Reader::readCollisions(size_t iteration) const {
  std::vector<HDF5Definitions::CollisionInfo> collisions{};
#ifdef LADDS_HDF5
  file.readDataset(
      collisions, HDF5Definitions::groupCollisionData + std::to_string(iteration) + HDF5Definitions::datasetCollisions);
#endif
  return collisions;
}

size_t HDF5Reader::readLastIterationNr() {
  std::vector<size_t> allIterations;
#ifdef LADDS_HDF5
  const auto groupStrLength = std::string(HDF5Definitions::groupParticleData).size();
  // gather all ID datasets. They contain the iteration number in their path
  auto allIdDatasets = file.findDatasets(HDF5Definitions::datasetParticleIDs);
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
