/**
 * @file VTUWriter.h
 * @author F. Gratl
 * @date 06.12.21
 */

#pragma once

#include "ladds/TypeDefinitions.h"
namespace VTUWriter {

/**
 * Write a vtk file with the current state of the simulation.
 */
void writeVTK(size_t iteration, const AutoPas_t &autopas);
void writeLegacyVTKBinary(size_t iteration, const AutoPas_t &autopas);

}  // namespace VTUWriter
