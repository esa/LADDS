/**
 * @file VTUWriter.h
 * @author F. Gratl
 * @date 06.12.21
 */

#include "VTUWriter.h"

#include <breakupModel/output/VTKWriter.h>

#include "ladds/particle/SatelliteToParticleConverter.h"

template <class T>
void swapEndianess(T &var) {
  char *varArray = reinterpret_cast<char *>(&var);
  for (long i = 0; i < static_cast<long>(sizeof(var) / 2); i++) {
    std::swap(varArray[sizeof(var) - 1 - i], varArray[i]);
  }
}

template <class T, bool swap, class F>
void printScalar(const AutoPas_t &autopas, std::ofstream &file, F fun) {
  for (const auto &particle : autopas) {
    T entry = fun(particle);
    if (swap) {
      swapEndianess(entry);
    }
    file.write(reinterpret_cast<char *>(&entry), sizeof(T));
  }
}

template <class T, bool swap, class F>
void printArray(const AutoPas_t &autopas, std::ofstream &file, F fun) {
  for (const auto &particle : autopas) {
    for (T entry : fun(particle)) {
      if (swap) {
        swapEndianess(entry);
      }
      file.write(reinterpret_cast<char *>(&entry), sizeof(T));
    }
  }
}

void VTUWriter::writeLegacyVTKBinary(size_t iteration, const AutoPas_t &autopas) {
  auto numParticles = autopas.getNumberOfParticles();

  std::ofstream file;
  file.open("test_" + std::to_string(iteration) + ".vtk", std::ios::out);

  file << "# vtk DataFile Version 2.0\n"
          "Timestep\n"
          "BINARY\n"
          "DATASET STRUCTURED_GRID\n"
          "DIMENSIONS 1 1 1\n";

  file << "POINTS " << numParticles << " double\n";
  printArray<float, true>(autopas, file, [](const auto &p) { return p.getR(); });
  file << "POINT_DATA " << numParticles
       << "\n"
          "VECTORS velocities double\n";
  printArray<float, true>(autopas, file, [](const auto &p) { return p.getV(); });
  file << "SCALARS particleIds Int32\n"
       << "LOOKUP_TABLE default\n";
  printScalar<unsigned short, true>(autopas, file, [](const auto &p) { return p.getID(); });

  file.close();
}

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
