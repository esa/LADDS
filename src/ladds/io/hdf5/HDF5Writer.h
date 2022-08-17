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

#include "HDF5Definitions.h"
#include "ladds/TypeDefinitions.h"
#include "ladds/io/ConjunctionWriterInterface.h"

namespace LADDS {

/**
 * Wrapper for the whole logic of writing the HDF5 file.
 * All data is written in 32 bit precision.
 */
class HDF5Writer final : public ConjuctionWriterInterface {
 public:
  /**
   * Constructor setting up the file and creating the custom data type.
   * @note Here the actual file is created on disk.
   * @param filename
   * @param replace if true replace an existing file, else append.
   * @param compressionLevel
   */
  HDF5Writer(const std::string &filename, bool replace, unsigned int compressionLevel);

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
   * Write evasions of the current iteration to the HDF5 file. If evasions is empty, nothing is written.
   * @param iteration
   * @param evasions
   */
  void writeEvasions(size_t iteration, const CollisionFunctor::CollisionCollectionT &evasions) override;

 private:
  /**
   * Highest partilce ID that was written in any previous iteration.
   * For anything below this ID constant particle properties are already written.
   */
  HDF5Definitions::IntType maxWrittenParticleID{0};
#ifdef LADDS_HDF5
  /**
   * Actual file that will be created. All of the data this writer gets ends up in this one file.
   */
  h5pp::File _file;
  /**
   * Object holding the info for the hdf5 compound type of the collision data.
   */
  h5pp::hid::h5t collisionInfoH5Type;
  /**
   * Object holding the info for the hdf5 compound type of the collision data.
   */
  h5pp::hid::h5t evasionInfoH5Type;
  /**
   * Object holding the info for the hdf5 compound type of the constant particle properties data.
   */
  h5pp::hid::h5t particleConstantPropertiesH5Type;
#endif
};
}  // namespace LADDS