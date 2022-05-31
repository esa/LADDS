/**
 * @file VTUWriter.h
 * @author F. Gratl
 * @date 06.12.21
 */

#include "VTUWriter.h"

#include <breakupModel/output/VTKWriter.h>

#include "ladds/particle/SatelliteToParticleConverter.h"

namespace LADDS {

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

void VTUWriter::writeVTU(ConfigReader &config,
                         size_t iteration,
                         const AutoPas_t &autopas,
                         const DomainDecomposition &decomposition) {
  VTKWriter vtkWriter(filenamePayload(config, iteration, decomposition));
  std::vector<Satellite> allParticles;
  allParticles.reserve(autopas.getNumberOfParticles());
  for (const auto &p : autopas) {
    auto sat = SatelliteToParticleConverter::convertParticleToSatellite(p);
    sat.setPosition(autopas::utils::ArrayMath::mulScalar(sat.getPosition(), 1 / 1000.));
    allParticles.push_back(sat);
  }
  // sort particles by Id to provide consistent output files
  std::sort(
      allParticles.begin(), allParticles.end(), [](const auto &p1, const auto &p2) { return p1.getId() < p2.getId(); });
  vtkWriter.printResult(allParticles);
}

void VTUWriter::writePVTU(ConfigReader &config, size_t iteration, const DomainDecomposition &decomposition) {
  const auto filename = filenameMetadata(config, iteration);
  std::ofstream pvtuFile;
  pvtuFile.open(filename, std::ios::out | std::ios::binary);

  if (not pvtuFile.is_open()) {
    throw std::runtime_error("Simulation::writeVTKFile(): Failed to open file \"" + filename + "\"");
  }
  pvtuFile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n";
  pvtuFile << "<VTKFile byte_order=\"LittleEndian\" type=\"PUnstructuredGrid\" version=\"0.1\">\n";
  pvtuFile << "  <PUnstructuredGrid>\n";
  pvtuFile << "    <PPointData>\n";
  pvtuFile << "      <PDataArray type=\"Float32\" Name=\"characteristic-length\" />\n";
  pvtuFile << "      <PDataArray type=\"Float32\" Name=\"mass\" />\n";
  pvtuFile << "      <PDataArray type=\"Float32\" Name=\"area\" />\n";
  pvtuFile << "      <PDataArray type=\"Float32\" Name=\"area-to-mass\" />\n";
  pvtuFile << "      <PDataArray type=\"Float32\" Name=\"velocity\" />\n";
  pvtuFile << "      <PDataArray type=\"Float32\" Name=\"ejection-velocity\" />\n";
  pvtuFile << "    </PPointData>\n";
  pvtuFile << "    <PCellData/>\n";
  pvtuFile << "    <PPoints>\n";
  pvtuFile << "      <PDataArray Name=\"position\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\" />\n";
  pvtuFile << "    </PPoints>\n";
  pvtuFile << "    <PCells>\n";
  pvtuFile << "      <PDataArray Name=\"types\" NumberOfComponents=\"0\" format=\"ascii\" type=\"Float32\"/>\n";
  pvtuFile << "    </PCells>\n";

  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(decomposition.getCommunicator(), &numRanks);
  for (int rank = 0; rank < numRanks; ++rank) {
    pvtuFile << "    <Piece Source=\"" << filenamePayload(config, iteration, decomposition, rank) << "\"/>\n";
  }

  pvtuFile << "  </PUnstructuredGrid>\n";
  pvtuFile << "</VTKFile>\n";

  pvtuFile.close();
}

std::string VTUWriter::filenameMetadata(ConfigReader &config, size_t iteration) {
  const auto maxDigitsIterations = std::to_string(config.getLastIterationNr()).size();
  std::stringstream ss;
  ss << "Output_" << std::setfill('0') << std::setw(static_cast<int>(maxDigitsIterations)) << iteration << ".pvtu";
  return ss.str();
}

std::string VTUWriter::filenamePayload(ConfigReader &config,
                                       size_t iteration,
                                       const DomainDecomposition &decomposition,
                                       int rank) {
  if (rank == -1) {
    autopas::AutoPas_MPI_Comm_rank(decomposition.getCommunicator(), &rank);
  }
  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(decomposition.getCommunicator(), &numRanks);
  const auto maxDigitsIterations = std::to_string(config.getLastIterationNr()).size();
  std::stringstream ss;
  ss << "Output_" << std::setfill('0') << std::setw(static_cast<int>(maxDigitsIterations)) << iteration << "_"
     << std::setw(static_cast<int>(std::to_string(numRanks).size())) << rank << ".vtu";
  return ss.str();
}

}  // namespace LADDS