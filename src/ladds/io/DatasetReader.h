#pragma once

#include <breakupModel/input/CSVReader.h>

#include "ladds/particle/Particle.h"

namespace DatasetReader {

/**
 * Reads the passed position and velocity csv files. Returns a vector of particles.
 */
std::vector<Particle> readDataset(const std::string &csvFilepath);
}  // namespace DatasetReader