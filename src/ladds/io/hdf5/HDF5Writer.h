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
