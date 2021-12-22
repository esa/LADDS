/**
 * @file HDF5Writer.h
 * @author F. Gratl
 * @date 14.12.21
 */

#pragma once

#include <h5pp/h5pp.h>

#include <string>

#include "ConjunctionWriterInterface.h"
#include "ladds/TypeDefinitions.h"

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
  void writeConjunctions(size_t iteration,
                         const std::unordered_map<Particle *, std::tuple<Particle *, double>> &collisions) override;

 private:
  /**
   * Actual file that will be created. All of the data this writer gets ends up in this one file.
   */
  h5pp::File _file;

  const std::string groupParticleData = "ParticleData/";

  const std::string groupCollisionData = "CollisionData/";

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
  };

  /*
   *
   */
  h5pp::hid::h5t collisionInfoH5Type;
};
