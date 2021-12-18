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
      : _file(filename, h5pp::FilePermission::REPLACE) {
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
  template <class T>
  struct Vec3 {
    T x, y, z;
  };

  /**
   * Actual file that will be created. All of the data this writer gets ends up in this one file.
   */
  h5pp::File _file;
};
