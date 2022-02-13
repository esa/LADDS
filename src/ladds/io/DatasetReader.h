#pragma once

#include <breakupModel/input/CSVReader.h>

#include "ladds/particle/Particle.h"

namespace DatasetReader {

/**
 * Reads the passed position and velocity csv files. Returns a vector of particles.
 * @param csvFilepath
 * @param coefficientOfDrag c_D used to initialize all loaded particles for which no bstar is available.
 * @return
 */
std::vector<Particle> readDataset(const std::string &csvFilepath, double coefficientOfDrag);
}  // namespace DatasetReader
