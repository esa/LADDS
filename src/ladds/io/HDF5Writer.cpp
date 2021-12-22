/**
 * @file HDF5Writer.cpp
 * @author F. Gratl
 * @date 14.12.21
 */

#include "HDF5Writer.h"

HDF5Writer::HDF5Writer(const std::string &filename, unsigned int compressionLevel)
    : _file(filename, h5pp::FilePermission::REPLACE),
      collisionInfoH5Type(H5Tcreate(H5T_COMPOUND, sizeof(CollisionInfo))) {
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
}

void HDF5Writer::writeParticles(size_t iteration, const AutoPas_t &autopas) {
  // Data will be casted to these types before writing
  using floatType = float;
  using intType = unsigned int;

  std::vector<Vec3<floatType>> vecPos;
  std::vector<Vec3<floatType>> vecVel;
  std::vector<intType> vecId;

  vecPos.reserve(autopas.getNumberOfParticles());
  vecVel.reserve(autopas.getNumberOfParticles());
  vecId.reserve(autopas.getNumberOfParticles());

  for (const auto &particle : autopas) {
    const auto &pos = particle.getR();
    const auto &vel = particle.getV();
    // pack data and make sure it is of the correct type
    vecPos.emplace_back<Vec3<floatType>>(
        {static_cast<floatType>(pos[0]), static_cast<floatType>(pos[1]), static_cast<floatType>(pos[2])});
    vecVel.emplace_back<Vec3<floatType>>(
        {static_cast<floatType>(vel[0]), static_cast<floatType>(vel[1]), static_cast<floatType>(vel[2])});
    vecId.emplace_back(static_cast<intType>(particle.getID()));
  }

  const auto group = groupParticleData + std::to_string(iteration);
  const auto &compressionLvl = _file.getCompressionLevel();
  _file.writeDataset_compressed(vecPos, group + "/Particles/Positions", compressionLvl);
  _file.writeDataset_compressed(vecVel, group + "/Particles/Velocities", compressionLvl);
  _file.writeDataset_compressed(vecId, group + "/Particles/IDs", compressionLvl);
}

void HDF5Writer::writeConjunctions(size_t iteration,
                                   const std::unordered_map<Particle *, std::tuple<Particle *, double>> &collisions) {
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
  _file.writeDataset(data, group + "/Collisions", collisionInfoH5Type);
}
