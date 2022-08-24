/**
 * @file HDF5Definitions.h
 * @author F. Gratl
 * @date 19.01.22
 */

#pragma once

#include <string>

namespace LADDS::HDF5Definitions {
constexpr auto groupParticleData = "ParticleData/";
constexpr auto datasetParticlePositions = "/Particles/Positions";
constexpr auto datasetParticleVelocities = "/Particles/Velocities";
constexpr auto datasetParticleIDs = "/Particles/IDs";
constexpr auto tableParticleConstantProperties = "ConstantProperties";
constexpr auto groupCollisionData = "CollisionData/";
constexpr auto datasetCollisions = "/Collisions";
constexpr auto groupEvasionData = "EvasionData/";
constexpr auto datasetEvasions = "/Evasions";

/**
 * Type to which any floating point data will be cast before writing.
 */
using FloatType = float;
/**
 * Type to which any integer data will be cast before writing.
 */
using IntType = unsigned int;

/**
 * This represents one line of 3D vector data in the HDF5 file.
 */
template <class T>
struct Vec3 {
  T x, y, z;
};

/**
 * Type for the information of a single collision.
 */
struct CollisionInfo {
  unsigned int idA, idB;
  float distanceSquared;
  bool operator==(const CollisionInfo &rhs) const {
    return idA == rhs.idA && idB == rhs.idB && distanceSquared == rhs.distanceSquared;
  }
  bool operator!=(const CollisionInfo &rhs) const {
    return !(rhs == *this);
  }
};

/**
 * Type for the information of a single evaded collision.
 */
struct EvasionInfo {
  unsigned int idA, idB;
  float distanceSquared;
  bool operator==(const EvasionInfo &rhs) const {
    return idA == rhs.idA && idB == rhs.idB && distanceSquared == rhs.distanceSquared;
  }
  bool operator!=(const EvasionInfo &rhs) const {
    return !(rhs == *this);
  }
};

/**
 * Type for information of a single particle that will stay constant throughout the simulation.
 */
struct ParticleConstantProperties {
  IntType id{};
  const char *identifier{};
  FloatType mass{}, radius{}, bcInv{};
  int activityState{};
};
}  // namespace LADDS::HDF5Definitions
