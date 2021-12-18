/**
 * @file HDF5Writer.cpp
 * @author F. Gratl
 * @date 14.12.21
 */

#include "HDF5Writer.h"

void HDF5Writer::write(size_t iteration, const AutoPas_t &autopas) {
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
        {static_cast<floatType>(pos[0]), static_cast<floatType>(pos[0]), static_cast<floatType>(pos[0])});
    vecVel.emplace_back<Vec3<floatType>>(
        {static_cast<floatType>(vel[0]), static_cast<floatType>(vel[0]), static_cast<floatType>(vel[0])});
    vecId.emplace_back(static_cast<intType>(particle.getID()));
  }

  const auto group = groupParticleData + std::to_string(iteration);
  const auto &compressionLvl = _file.getCompressionLevel();
  _file.writeDataset_compressed(vecPos, group + "/Particles/Positions", compressionLvl);
  _file.writeDataset_compressed(vecVel, group + "/Particles/Velocities", compressionLvl);
  _file.writeDataset_compressed(vecId, group + "/Particles/IDs", compressionLvl);
}
