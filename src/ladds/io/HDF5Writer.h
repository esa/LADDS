/**
 * @file HDF5Writer.h
 * @author F. Gratl
 * @date 14.12.21
 */

#pragma once

#include <h5pp/h5pp.h>

#include <string>

#include "ladds/TypeDefinitions.h"

class HDF5Writer {
 public:
  explicit HDF5Writer(const std::string &filename) : _file(filename, h5pp::FilePermission::REPLACE) {}

  void write(size_t iteration, const AutoPas_t &autopas);

 private:
  h5pp::File _file;
};
