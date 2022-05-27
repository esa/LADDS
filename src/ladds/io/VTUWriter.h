/**
 * @file VTUWriter.h
 * @author F. Gratl
 * @date 06.12.21
 */

#pragma once

#include "ladds/TypeDefinitions.h"
#include "ladds/distributedMemParallelization/DomainDecomposition.h"
#include "ladds/io/ConfigReader.h"

namespace LADDS::VTUWriter {

/**
 * Write a vtu file (xml style) with the current state of the simulation using the VTKWriter from the breakup code.
 *
 * @param config
 * @param iteration
 * @param autopas
 * @param decomposition
 */
void writeVTU(ConfigReader &config,
              size_t iteration,
              const AutoPas_t &autopas,
              const DomainDecomposition &decomposition);

/**
 * Write a pvtu file which references all vtu files from all ranks in time steps
 * @param iteration
 * @param decomposition
 */
void writePVTU(ConfigReader &config, size_t iteration, const DomainDecomposition &decomposition);

/**
 * Write a vtk file in legacy table style (no xml) but with binary storage layout.
 * Without compression this is probably the smallest output file possible.
 * @param iteration
 * @param autopas
 */
void writeLegacyVTKBinary(size_t iteration, const AutoPas_t &autopas);

/**
 * Generates a filename for the metadata file (pvtu).
 * @param iteration
 * @return
 */
std::string filenameMetadata(ConfigReader &config, size_t iteration);

/**
 * Generates a filename for the payload file (vtu) of a specific rank.
 * @param iteration
 * @param rank If omitted, the current rank number is chosen.
 * @return
 */
std::string filenamePayload(ConfigReader &config,
                            size_t iteration,
                            const DomainDecomposition &decomposition,
                            int rank = -1);
}  // namespace LADDS::VTUWriter
