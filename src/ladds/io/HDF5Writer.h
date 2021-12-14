/**
 * @file HDF5Writer.h
 * @author F. Gratl
 * @date 14.12.21
 */

#pragma once

#include <highfive/H5File.hpp>
#include <string>

#include "ladds/TypeDefinitions.h"

class HDF5Writer {
 public:
  explicit HDF5Writer(const std::string &filename) : _file(filename, HighFive::File::Overwrite) {
    //    HighFive::CompoundType compoundType{{"x", HighFive::AtomicType<double>{}},
    //                                        {"y", HighFive::AtomicType<double>{}},
    //                                        {"z", HighFive::AtomicType<double>{}}};
    //    compoundType.commit(_file, "PositionsType");

    //    HighFive::CompoundType array3D_t({{"x", HighFive::AtomicType<double>{}},
    //                                      {"y", HighFive::AtomicType<double>{}},
    //                                      {"z", HighFive::AtomicType<double>{}}});
    //    HighFive::CompoundType positionsVector_t(std::vector<HighFive::CompoundType::member_def>{array3D_t{}});
    //    positionsVector_t.commit(_file, "PositionsType");
  }

  void write(size_t iteration, const AutoPas_t &autopas);

 private:
  struct data_t {
    size_t id;
    std::array<double, 3> r;
    std::array<double, 3> v;
  } __attribute__((aligned(64)));

  HighFive::File _file;
};
