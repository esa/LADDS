/**
 * @file VTUWriter.h
 * @author F. Gratl
 * @date 06.12.21
 */

#include "VTUWriter.h"

#include <breakupModel/output/VTKWriter.h>

#include "ladds/particle/SatelliteToParticleConverter.h"

void VTUWriter::writeVTK(size_t iteration, const AutoPas_t &autopas) {
  VTKWriter vtkWriter("output_" + std::to_string(iteration) + ".vtu");
  std::vector<Satellite> allParticles;
  allParticles.reserve(autopas.getNumberOfParticles());
  for (const auto &p : autopas) {
    allParticles.push_back(SatelliteToParticleConverter::convertParticleToSatellite(p));
  }
  // sort particles by Id to provide consistent output files
  std::sort(
      allParticles.begin(), allParticles.end(), [](const auto &p1, const auto &p2) { return p1.getId() < p2.getId(); });
  vtkWriter.printResult(allParticles);
}
