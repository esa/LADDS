/**
 * @file HDF5Reader.h
 * @author F. Gratl
 * @date 23.12.21
 */

#pragma once

#ifdef LADDS_HDF5
#include <h5pp/h5pp.h>
#endif

#include <string>

#include "ladds/particle/Particle.h"
class HDF5Reader {
 public:
  explicit HDF5Reader(const std::string &filename);

 private:
  std::vector<Particle> readParticles(size_t iteration);

#ifdef LADDS_HDF5
  /**
   * Actual file that is read.
   */
  h5pp::File file;
#endif
};
