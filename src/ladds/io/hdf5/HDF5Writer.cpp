/**
 * @file HDF5Writer.cpp
 * @author F. Gratl
 * @date 14.12.21
 */

#include "HDF5Writer.h"

#include "HDF5Definitions.h"

HDF5Writer::HDF5Writer(const std::string &filename, unsigned int compressionLevel)
#ifdef LADDS_HDF5
    : _file(filename, h5pp::FilePermission::REPLACE),
      collisionInfoH5Type(H5Tcreate(H5T_COMPOUND, sizeof(HDF5Definitions::CollisionInfo)))
#endif
{
#ifdef LADDS_HDF5
  _file.setCompressionLevel(compressionLevel);
  H5Tinsert(collisionInfoH5Type,
            "idA",
            HOFFSET(HDF5Definitions::CollisionInfo, idA),
            h5pp::type::getH5NativeType<decltype(HDF5Definitions::CollisionInfo::idA)>());
  H5Tinsert(collisionInfoH5Type,
            "idB",
            HOFFSET(HDF5Definitions::CollisionInfo, idB),
            h5pp::type::getH5NativeType<decltype(HDF5Definitions::CollisionInfo::idB)>());
  H5Tinsert(collisionInfoH5Type,
            "distanceSquared",
            HOFFSET(HDF5Definitions::CollisionInfo, distanceSquared),
            h5pp::type::getH5NativeType<decltype(HDF5Definitions::CollisionInfo::distanceSquared)>());
#else
  throw std::runtime_error("LADDS was compiled without HDF5 support, so the HDF5Writer can't do anything!");
#endif
}

void HDF5Writer::writeParticles(size_t iteration, const AutoPas_t &autopas) {
#ifdef LADDS_HDF5

  std::vector<HDF5Definitions::Vec3<HDF5Definitions::FloatType>> vecPos;
  std::vector<HDF5Definitions::Vec3<HDF5Definitions::FloatType>> vecVel;
  std::vector<HDF5Definitions::IntType> vecId;

  vecPos.reserve(autopas.getNumberOfParticles());
  vecVel.reserve(autopas.getNumberOfParticles());
  vecId.reserve(autopas.getNumberOfParticles());

  for (const auto &particle : autopas) {
    const auto &pos = particle.getR();
    const auto &vel = particle.getV();
    // pack data and make sure it is of the correct type
    vecPos.emplace_back<HDF5Definitions::Vec3<HDF5Definitions::FloatType>>(
        {static_cast<HDF5Definitions::FloatType>(pos[0]),
         static_cast<HDF5Definitions::FloatType>(pos[1]),
         static_cast<HDF5Definitions::FloatType>(pos[2])});
    vecVel.emplace_back<HDF5Definitions::Vec3<HDF5Definitions::FloatType>>(
        {static_cast<HDF5Definitions::FloatType>(vel[0]),
         static_cast<HDF5Definitions::FloatType>(vel[1]),
         static_cast<HDF5Definitions::FloatType>(vel[2])});
    vecId.emplace_back(static_cast<HDF5Definitions::IntType>(particle.getID()));
  }

  const auto group = HDF5Definitions::groupParticleData + std::to_string(iteration);
  const auto &compressionLvl = _file.getCompressionLevel();
  _file.writeDataset_compressed(vecPos, group + HDF5Definitions::datasetParticlePositions, compressionLvl);
  _file.writeDataset_compressed(vecVel, group + HDF5Definitions::datasetParticleVelocities, compressionLvl);
  _file.writeDataset_compressed(vecId, group + HDF5Definitions::datasetParticleIDs, compressionLvl);
#endif
}

void HDF5Writer::writeConjunctions(size_t iteration, const CollisionFunctor::CollisionCollectionT &collisions) {
#ifdef LADDS_HDF5
  if (collisions.empty()) {
    return;
  }
  std::vector<HDF5Definitions::CollisionInfo> data;
  data.reserve(collisions.size());

  for (const auto &[p1, p2, distanceSquare] : collisions) {
    data.emplace_back<HDF5Definitions::CollisionInfo>(
        {static_cast<decltype(HDF5Definitions::CollisionInfo::idA)>(p1->getID()),
         static_cast<decltype(HDF5Definitions::CollisionInfo::idB)>(p2->getID()),
         static_cast<decltype(HDF5Definitions::CollisionInfo::distanceSquared)>(distanceSquare)});
  }
  const auto group = HDF5Definitions::groupCollisionData + std::to_string(iteration);
  _file.writeDataset(data, group + HDF5Definitions::datasetCollisions, collisionInfoH5Type);
#endif
}
