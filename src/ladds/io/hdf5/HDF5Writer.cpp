/**
 * @file HDF5Writer.cpp
 * @author F. Gratl
 * @date 14.12.21
 */

#include "HDF5Writer.h"

HDF5Writer::HDF5Writer(const std::string &filename, unsigned int compressionLevel)
#ifdef LADDS_HDF5
    : _file(filename, h5pp::FilePermission::REPLACE),
      collisionInfoH5Type(H5Tcreate(H5T_COMPOUND, sizeof(CollisionInfo)))
#endif
{
#ifdef LADDS_HDF5
  _file.setCompressionLevel(compressionLevel);
  H5Tinsert(collisionInfoH5Type,
            "idA",
            HOFFSET(CollisionInfo, idA),
            h5pp::type::getH5NativeType<decltype(CollisionInfo::idA)>());
  H5Tinsert(collisionInfoH5Type,
            "idB",
            HOFFSET(CollisionInfo, idB),
            h5pp::type::getH5NativeType<decltype(CollisionInfo::idB)>());
  H5Tinsert(collisionInfoH5Type,
            "distanceSquared",
            HOFFSET(CollisionInfo, distanceSquared),
            h5pp::type::getH5NativeType<decltype(CollisionInfo::distanceSquared)>());
#else
  throw std::runtime_error("LADDS was compiled without HDF5 support, so the HDF5Writer can't do anything!");
#endif
}

void HDF5Writer::writeParticles(size_t iteration, const AutoPas_t &autopas) {
#ifdef LADDS_HDF5

  std::vector<Vec3<FloatType>> vecPos;
  std::vector<Vec3<FloatType>> vecVel;
  std::vector<IntType> vecId;

  vecPos.reserve(autopas.getNumberOfParticles());
  vecVel.reserve(autopas.getNumberOfParticles());
  vecId.reserve(autopas.getNumberOfParticles());

  for (const auto &particle : autopas) {
    const auto &pos = particle.getR();
    const auto &vel = particle.getV();
    // pack data and make sure it is of the correct type
    vecPos.emplace_back<Vec3<FloatType>>(
        {static_cast<FloatType>(pos[0]), static_cast<FloatType>(pos[1]), static_cast<FloatType>(pos[2])});
    vecVel.emplace_back<Vec3<FloatType>>(
        {static_cast<FloatType>(vel[0]), static_cast<FloatType>(vel[1]), static_cast<FloatType>(vel[2])});
    vecId.emplace_back(static_cast<IntType>(particle.getID()));
  }

  const auto group = groupParticleData + std::to_string(iteration);
  const auto &compressionLvl = _file.getCompressionLevel();
  _file.writeDataset_compressed(vecPos, group + datasetParticlePositions, compressionLvl);
  _file.writeDataset_compressed(vecVel, group + datasetParticleVelocities, compressionLvl);
  _file.writeDataset_compressed(vecId, group + datasetParticleIDs, compressionLvl);
#endif
}

void HDF5Writer::writeConjunctions(size_t iteration,
                                   const std::unordered_map<Particle *, std::tuple<Particle *, double>> &collisions) {
#ifdef LADDS_HDF5
  if (collisions.empty()) {
    return;
  }
  std::vector<CollisionInfo> data;
  data.reserve(collisions.size());

  for (const auto &[p1, p2AndDistanceSquare] : collisions) {
    const auto &[p2, distanceSquare] = p2AndDistanceSquare;
    data.emplace_back<CollisionInfo>({static_cast<decltype(CollisionInfo::idA)>(p1->getID()),
                                      static_cast<decltype(CollisionInfo::idB)>(p2->getID()),
                                      static_cast<decltype(CollisionInfo::distanceSquared)>(distanceSquare)});
  }
  const auto group = groupCollisionData + std::to_string(iteration);
  _file.writeDataset(data, group + datasetCollisions, collisionInfoH5Type);
#endif
}

bool HDF5Writer::CollisionInfo::operator==(const HDF5Writer::CollisionInfo &rhs) const {
  return idA == rhs.idA && idB == rhs.idB && distanceSquared == rhs.distanceSquared;
}
bool HDF5Writer::CollisionInfo::operator!=(const HDF5Writer::CollisionInfo &rhs) const {
  return !(rhs == *this);
}
