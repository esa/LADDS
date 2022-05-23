/**
 * @file HDF5Writer.cpp
 * @author F. Gratl
 * @date 14.12.21
 */

#include "HDF5Writer.h"

#include "HDF5Definitions.h"

HDF5Writer::HDF5Writer(const std::string &filename, bool replace, unsigned int compressionLevel)
#ifdef LADDS_HDF5
    : _file(filename, replace ? h5pp::FilePermission::REPLACE : h5pp::FilePermission::READWRITE),
      collisionInfoH5Type(H5Tcreate(H5T_COMPOUND, sizeof(HDF5Definitions::CollisionInfo))),
      particleConstantPropertiesH5Type(H5Tcreate(H5T_COMPOUND, sizeof(HDF5Definitions::ParticleConstantProperties)))
#endif
{
#ifdef LADDS_HDF5
  // if the file already exists we don't need to set a compression level
  if (replace) {
    _file.setCompressionLevel(compressionLevel);
  }
  // CollisionInfo
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

  // ParticleConstantProperties
  H5Tinsert(particleConstantPropertiesH5Type,
            "id",
            HOFFSET(HDF5Definitions::ParticleConstantProperties, id),
            h5pp::type::getH5NativeType<decltype(HDF5Definitions::ParticleConstantProperties::id)>());
  // Define variable length type for c-string field
  h5pp::hid::h5t H5_VLEN_STR_TYPE = H5Tcopy(H5T_C_S1);
  H5Tset_size(H5_VLEN_STR_TYPE, H5T_VARIABLE);
  H5Tset_strpad(H5_VLEN_STR_TYPE, H5T_STR_NULLTERM);
  H5Tset_cset(H5_VLEN_STR_TYPE, H5T_CSET_UTF8);
  H5Tinsert(particleConstantPropertiesH5Type,
            "identifier",
            HOFFSET(HDF5Definitions::ParticleConstantProperties, identifier),
            H5_VLEN_STR_TYPE);
  H5Tinsert(particleConstantPropertiesH5Type,
            "mass",
            HOFFSET(HDF5Definitions::ParticleConstantProperties, mass),
            h5pp::type::getH5NativeType<decltype(HDF5Definitions::ParticleConstantProperties::mass)>());
  H5Tinsert(particleConstantPropertiesH5Type,
            "radius",
            HOFFSET(HDF5Definitions::ParticleConstantProperties, radius),
            h5pp::type::getH5NativeType<decltype(HDF5Definitions::ParticleConstantProperties::radius)>());
  H5Tinsert(particleConstantPropertiesH5Type,
            "bcInv",
            HOFFSET(HDF5Definitions::ParticleConstantProperties, bcInv),
            h5pp::type::getH5NativeType<decltype(HDF5Definitions::ParticleConstantProperties::bcInv)>());
  H5Tinsert(particleConstantPropertiesH5Type,
            "activityState",
            HOFFSET(HDF5Definitions::ParticleConstantProperties, activityState),
            h5pp::type::getH5NativeType<decltype(HDF5Definitions::ParticleConstantProperties::activityState)>());
#else
  throw std::runtime_error("LADDS was compiled without HDF5 support, so the HDF5Writer can't do anything!");
#endif
}

void HDF5Writer::writeParticles(size_t iteration, const AutoPas_t &autopas) {
#ifdef LADDS_HDF5

  std::vector<HDF5Definitions::Vec3<HDF5Definitions::FloatType>> vecPos;
  std::vector<HDF5Definitions::Vec3<HDF5Definitions::FloatType>> vecVel;
  std::vector<HDF5Definitions::IntType> vecId;
  HDF5Definitions::IntType maxParticleId{0};
  std::vector<HDF5Definitions::ParticleConstantProperties> newConstantProperties;

  vecPos.reserve(autopas.getNumberOfParticles());
  vecVel.reserve(autopas.getNumberOfParticles());
  vecId.reserve(autopas.getNumberOfParticles());

  for (const auto &particle : autopas) {
    const auto &pos = particle.getR();
    const auto &vel = particle.getV();
    const auto &id = static_cast<HDF5Definitions::IntType>(particle.getID());
    // pack data and make sure it is of the correct type
    vecPos.emplace_back<HDF5Definitions::Vec3<HDF5Definitions::FloatType>>(
        {static_cast<HDF5Definitions::FloatType>(pos[0]),
         static_cast<HDF5Definitions::FloatType>(pos[1]),
         static_cast<HDF5Definitions::FloatType>(pos[2])});
    vecVel.emplace_back<HDF5Definitions::Vec3<HDF5Definitions::FloatType>>(
        {static_cast<HDF5Definitions::FloatType>(vel[0]),
         static_cast<HDF5Definitions::FloatType>(vel[1]),
         static_cast<HDF5Definitions::FloatType>(vel[2])});
    vecId.emplace_back(id);
    // track the highest particle id that was written to the file
    // All particles that have a higher id than the highest id from the last time something was written are new
    // and their static properties need to be recorded.
    maxParticleId = std::max(maxParticleId, id);
    if (maxWrittenParticleID == 0 or id > maxWrittenParticleID) {
      newConstantProperties.emplace_back(
          HDF5Definitions::ParticleConstantProperties{id,
                                                      particle.getIdentifier().c_str(),
                                                      static_cast<HDF5Definitions::FloatType>(particle.getMass()),
                                                      static_cast<HDF5Definitions::FloatType>(particle.getRadius()),
                                                      static_cast<HDF5Definitions::FloatType>(particle.getBcInv()),
                                                      static_cast<int>(particle.getActivityState())});
    }
  }

  const auto group = HDF5Definitions::groupParticleData + std::to_string(iteration);
  const auto &compressionLvl = _file.getCompressionLevel();
  _file.writeDataset_compressed(vecPos, group + HDF5Definitions::datasetParticlePositions, compressionLvl);
  _file.writeDataset_compressed(vecVel, group + HDF5Definitions::datasetParticleVelocities, compressionLvl);
  _file.writeDataset_compressed(vecId, group + HDF5Definitions::datasetParticleIDs, compressionLvl);

  if (not newConstantProperties.empty()) {
    const auto particleConstantPropertiesFullPath =
        std::string(HDF5Definitions::groupParticleData) + HDF5Definitions::tableParticleConstantProperties;
    // it table does not exist yet create it
    if (not _file.linkExists(particleConstantPropertiesFullPath)) {
      _file.createTable(particleConstantPropertiesH5Type, particleConstantPropertiesFullPath, "ConstantProperties");
    }
    _file.appendTableRecords(newConstantProperties, particleConstantPropertiesFullPath);
  }

  maxWrittenParticleID = maxParticleId;
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
