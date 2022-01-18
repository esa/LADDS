/**
 * @file HDF5Writer.h
 * @author F. Gratl
 * @date 14.12.21
 */

#pragma once

#ifdef LADDS_HDF5
#include <h5pp/h5pp.h>
#endif

#include <string>

#include "ladds/TypeDefinitions.h"
#include "ladds/io/ConjunctionWriterInterface.h"

/**
 * Wrapper for the whole logic of writing the HDF5 file.
 * All data is written in 32 bit precision.
 */
class HDF5Writer final : public ConjuctionWriterInterface {
 public:
  /**
   * Constructor setting up the file and creating the custom data type.
   * @param filename
   * @param compressionLevel
   */
  HDF5Writer(const std::string &filename, unsigned int compressionLevel);

  ~HDF5Writer() override = default;

  /**
   * Write one dataset of all particle data in the current iteration.
   * @param iteration
   * @param autopas
   */
  void writeParticles(size_t iteration, const AutoPas_t &autopas);

  /**
   * Write collisions of the current iteration to the HDF5 file. If collisions is empty, nothing is written.
   * @param iteration
   * @param collisions
   */
  void writeConjunctions(size_t iteration, const CollisionFunctor::CollisionCollectionT &collisions) override;

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

  static constexpr auto groupParticleData = "ParticleData/";
  static constexpr auto datasetParticlePositions = "/Particles/Positions";
  static constexpr auto datasetParticleVelocities = "/Particles/Velocities";
  static constexpr auto datasetParticleIDs = "/Particles/IDs";

  static constexpr auto groupCollisionData = "CollisionData/";
  static constexpr auto datasetCollisions = "/Collisions";

  /**
   * Type for the information of a single collision.
   */
  struct CollisionInfo {
    unsigned int idA, idB;
    float distanceSquared;
    bool operator==(const CollisionInfo &rhs) const;
    bool operator!=(const CollisionInfo &rhs) const;
  };

 private:
#ifdef LADDS_HDF5
  /**
   * Actual file that will be created. All of the data this writer gets ends up in this one file.
   */
  h5pp::File _file;
#endif

#ifdef LADDS_HDF5
  /**
   * Object holding the info for the hdf5 compound type of the collision data.
   */
  h5pp::hid::h5t collisionInfoH5Type;
#endif
};
