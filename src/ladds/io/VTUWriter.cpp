/**
 * @file VTUWriter.h
 * @author F. Gratl
 * @date 06.12.21
 */

#include "VTUWriter.h"

#include <breakupModel/output/VTKWriter.h>

#include "ladds/particle/SatelliteToParticleConverter.h"

namespace LADDS {

VTUWriter::VTUWriter(ConfigReader &config, size_t iteration, const DomainDecomposition &decomposition)
    : _loggerName{VTUWriter::filenamePayload(config, iteration, decomposition)} {
  auto logger{spdlog::basic_logger_mt<spdlog::async_factory>(_loggerName, _loggerName)};
  logger->set_pattern("%v");
}
VTUWriter::~VTUWriter() {
  spdlog::drop(_loggerName);
}

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

void VTUWriter::writeVTU(const AutoPas_t &autopas) {
  // Header
  this->printHeader(autopas.getNumberOfParticles());

  // Point properties
  int myRank{};
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &myRank);
  this->printProperty<int, Particle>("rank", [&](const Particle &p) {return myRank;}, autopas);
  this->printProperty<unsigned int, Particle>("ID", &Particle::getID, autopas);
  this->printProperty<Particle::ActivityState, Particle>("activityState", &Particle::getActivityState, autopas);
  this->printProperty<double, Particle>("mass", &Particle::getMass, autopas);
  this->printProperty<double, Particle>("radius", &Particle::getRadius, autopas);
  //  this->printProperty<double, Particle>("areaToMass", &Particle::getAom, autopas);
  //  this->printProperty<double, Particle>("bcInv", &Particle::getBcInv, autopas);
  //  this->printProperty<double, Particle>("absSpeed", &Particle::getSpeed, autopas);
  //  this->printProperty<double, Particle>("heightAboveGround", &Particle::getHeight, autopas);
  this->printProperty<std::array<double, 3>, Particle>("velocity", &Particle::getVelocity, autopas);
  //  this->printProperty<std::array<double, 3>, Particle>("accT0", &Particle::getAccT0, autopas);
  //  this->printProperty<std::array<double, 3>, Particle>("accT1", &Particle::getAccT1, autopas);

  // Separator between point and point-tore-cell data
  this->printSeparator();

  // Point properties related to cell (position)
  this->printProperty<std::array<double, 3>, Particle>("position", &Particle::getPosition, autopas);

  // Footer
  this->printFooter();
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
  pvtuFile << "      <PDataArray type=\"Int32\" Name=\"rank\" />\n";
  pvtuFile << "      <PDataArray type=\"Int32\" Name=\"ID\" />\n";
  pvtuFile << "      <PDataArray type=\"Int32\" Name=\"activityState\" />\n";
  pvtuFile << "      <PDataArray type=\"Float32\" Name=\"mass\" />\n";
  pvtuFile << "      <PDataArray type=\"Float32\" Name=\"radius\" />\n";
  //  pvtuFile << "      <PDataArray type=\"Float32\" Name=\"areaToMass\" />\n";
  //  pvtuFile << "      <PDataArray type=\"Float32\" Name=\"bcInv\" />\n";
  //  pvtuFile << "      <PDataArray type=\"Float32\" Name=\"absSpeed\" />\n";
  //  pvtuFile << "      <PDataArray type=\"Float32\" Name=\"heightAboveGround\" />\n";
  pvtuFile << "      <PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\" />\n";
  //  pvtuFile << "      <PDataArray type=\"Float32\" Name=\"accT0\" NumberOfComponents=\"3\" format=\"ascii\" />\n";
  //  pvtuFile << "      <PDataArray type=\"Float32\" Name=\"accT1\" NumberOfComponents=\"3\" format=\"ascii\" />\n";
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
  ss << "Output_" << std::setfill('0') << std::setw(static_cast<int>(std::to_string(numRanks).size())) << rank << "_"
     << std::setw(static_cast<int>(maxDigitsIterations)) << iteration << ".vtu";
  return ss.str();
}

void VTUWriter::printHeader(size_t size) const {
  auto logger = spdlog::get(_loggerName);
  logger->info(R"(<?xml version="1.0" encoding="UTF-8" standalone="no" ?>)");
  logger->info(R"(<VTKFile byte_order="LittleEndian" type="UnstructuredGrid" version="0.1">)");
  logger->info(R"(  <UnstructuredGrid>)");
  logger->info(R"(    <Piece NumberOfCells="0" NumberOfPoints="{}">)", size);
  logger->info(R"(      <PointData>)");
}

void VTUWriter::printSeparator() const {
  auto logger = spdlog::get(_loggerName);
  logger->info(R"(      </PointData>)");
  logger->info(R"(      <CellData/>)");
  logger->info(R"(      <Points>)");
}

void VTUWriter::printFooter() const {
  auto logger = spdlog::get(_loggerName);
  logger->info(R"(      </Points>)");
  logger->info(R"(      <Cells>)");
  logger->info(R"(        <DataArray Name="types" NumberOfComponents="0" format="ascii" type="Float32"/>)");
  logger->info(R"(      </Cells>)");
  logger->info(R"(    </Piece>)");
  logger->info(R"(  </UnstructuredGrid>)");
  logger->info(R"(</VTKFile>)");
}
}  // namespace LADDS