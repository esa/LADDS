/**
 * @file HDF5Writer.h
 * @author F. Gratl
 * @date 14.12.21
 */

#pragma once

#include <h5pp/h5pp.h>

#include <string>

#include "ladds/TypeDefinitions.h"

/**
 * Wrapper for the whole logic of writing the HDF5 file.
 * All data is written in 32 bit precision.
 */
class HDF5Writer {
 public:
  /**
   * Constructor setting up the file and creating the custom data type.
   * @param filename
   * @param compressionLevel
   */
  HDF5Writer(const std::string &filename, unsigned int compressionLevel)
      : _file(filename, h5pp::FilePermission::REPLACE),
        ParticleDataH5Type(H5Tcreate(H5T_COMPOUND, sizeof(ParticleData))) {
    // Register ParticleData type with HDF5
    H5Tinsert(ParticleDataH5Type, "rx", HOFFSET(ParticleData, rx), H5T_NATIVE_FLOAT);
    H5Tinsert(ParticleDataH5Type, "ry", HOFFSET(ParticleData, ry), H5T_NATIVE_FLOAT);
    H5Tinsert(ParticleDataH5Type, "rz", HOFFSET(ParticleData, rz), H5T_NATIVE_FLOAT);
    H5Tinsert(ParticleDataH5Type, "vx", HOFFSET(ParticleData, vx), H5T_NATIVE_FLOAT);
    H5Tinsert(ParticleDataH5Type, "vy", HOFFSET(ParticleData, vy), H5T_NATIVE_FLOAT);
    H5Tinsert(ParticleDataH5Type, "vz", HOFFSET(ParticleData, vz), H5T_NATIVE_FLOAT);
    H5Tinsert(ParticleDataH5Type, "id", HOFFSET(ParticleData, id), H5T_NATIVE_UINT);
    _file.setCompressionLevel(compressionLevel);
  }

  /**
   * Write one dataset of all particle data in the current iteration.
   * @param iteration
   * @param autopas
   */
  void write(size_t iteration, const AutoPas_t &autopas);

 private:
  const std::string groupParticleData = "ParticleData/";

  /**
   * This represents one line of data in the HDF5 file.
   */
  struct ParticleData {
    float rx, ry, rz;
    float vx, vy, vz;
    unsigned int id;
  };

  /**
   * Actual file that will be created. All of the data this writer gets ends up in this one file.
   */
  h5pp::File _file;
  /**
   * Object representing the wrapper type HDF5 uses to understand the non-trivial struct ParticleData.
   */
  h5pp::hid::h5t ParticleDataH5Type;
};
