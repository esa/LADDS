/**
 * @file RegularGridDecompositionLogger.cpp
 * @author F. Gratl
 * @date 24.05.22
 */

#include "RegularGridDecompositionLogger.h"

#include <autopas/utils/WrapMPI.h>

#include <fstream>
#include <iomanip>

LADDS::RegularGridDecompositionLogger::RegularGridDecompositionLogger(
    ConfigReader &config, const LADDS::RegularGridDecomposition &decomposition)
    : loggerName("RegularGridDecompositionLogger"),
      decomposition(decomposition),
      maxDigitsIterations{std::to_string(config.getLastIterationNr()).size()} {}

void LADDS::RegularGridDecompositionLogger::writeMetafile(const size_t iteration) const {
  const auto filename = filenameMetadata(iteration);

  std::ofstream pvtsFile;
  pvtsFile.open(filename, std::ios::out);

  if (not pvtsFile.is_open()) {
    throw std::runtime_error("Simulation::writeVTKFile(): Failed to open file \"" + filename + "\"");
  }

  const auto &[dimensions, periods, _] = decomposition.getGridInfo();

  const std::array<double, 3> globalBoxMin = decomposition.getGlobalBoxMin();
  const std::array<double, 3> globalBoxMax = decomposition.getGlobalBoxMax();
  pvtsFile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n";
  pvtsFile << "<VTKFile byte_order=\"LittleEndian\" type=\"PStructuredGrid\" version=\"0.1\">\n";
  pvtsFile << "  <PStructuredGrid WholeExtent=\""
           << "0 " << dimensions[0] << " 0 " << dimensions[1] << " 0 " << dimensions[2] << "\" GhostLevel=\"0\">\n";
  pvtsFile << "    <PPointData/>\n";
  pvtsFile << "    <PCellData>\n";
  pvtsFile << "      <PDataArray type=\"Float32\" Name=\"CellSizeFactor\" />\n";
  pvtsFile << "      <PDataArray type=\"Int32\" Name=\"Container\" />\n";
  pvtsFile << "      <PDataArray type=\"Int32\" Name=\"DataLayout\" />\n";
  pvtsFile << "      <PDataArray type=\"Int32\" Name=\"LoadEstimator\" />\n";
  pvtsFile << "      <PDataArray type=\"Int32\" Name=\"Traversal\" />\n";
  pvtsFile << "      <PDataArray type=\"Int32\" Name=\"Newton3\" />\n";
  pvtsFile << "      <PDataArray type=\"Int32\" Name=\"Rank\" />\n";
  pvtsFile << "    </PCellData>\n";
  pvtsFile << "    <PPoints>\n";
  pvtsFile << "      <DataArray NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">\n";
  pvtsFile << "        " << globalBoxMin[0] << " " << globalBoxMin[1] << " " << globalBoxMin[2] << "\n";
  pvtsFile << "        " << globalBoxMax[0] << " " << globalBoxMin[1] << " " << globalBoxMin[2] << "\n";
  pvtsFile << "        " << globalBoxMin[0] << " " << globalBoxMax[1] << " " << globalBoxMin[2] << "\n";
  pvtsFile << "        " << globalBoxMax[0] << " " << globalBoxMax[1] << " " << globalBoxMin[2] << "\n";
  pvtsFile << "        " << globalBoxMin[0] << " " << globalBoxMin[1] << " " << globalBoxMax[2] << "\n";
  pvtsFile << "        " << globalBoxMax[0] << " " << globalBoxMin[1] << " " << globalBoxMax[2] << "\n";
  pvtsFile << "        " << globalBoxMin[0] << " " << globalBoxMax[1] << " " << globalBoxMax[2] << "\n";
  pvtsFile << "        " << globalBoxMax[0] << " " << globalBoxMax[1] << " " << globalBoxMax[2] << "\n";
  pvtsFile << "      </DataArray>\n";
  pvtsFile << "    </PPoints>\n";

  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(decomposition.getCommunicator(), &numRanks);
  for (int rank = 0; rank < numRanks; ++rank) {
    const auto filenamePiece = filenamePayload(iteration, rank);
    std::array<int, 3> rankGridCoords{};
    autopas::AutoPas_MPI_Cart_coords(decomposition.getCommunicator(), rank, dimensions.size(), rankGridCoords.data());
    pvtsFile << "    <Piece "
             << "Extent=\"" << (rankGridCoords[0]) << " " << (rankGridCoords[0] + 1) << " " << (rankGridCoords[1])
             << " " << (rankGridCoords[1] + 1) << " " << (rankGridCoords[2]) << " " << (rankGridCoords[2] + 1) << "\" "
             << "Source=\"" << filenamePiece << "\"/>\n";
  }

  pvtsFile << "  </PStructuredGrid>\n";
  pvtsFile << "</VTKFile>\n";

  pvtsFile.close();
}

void LADDS::RegularGridDecompositionLogger::writePayload(const size_t iteration,
                                                         const autopas::Configuration &autoPasConfig) const {
  const auto timestepFileName = filenamePayload(iteration);

  std::ofstream vtsFile;
  vtsFile.open(timestepFileName, std::ios::out);

  const auto &[dimensions, periods, rankGridCoords] = decomposition.getGridInfo();
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(decomposition.getCommunicator(), &rank);
  const std::array<double, 3> localBoxMin = decomposition.getLocalBoxMin();
  const std::array<double, 3> localBoxMax = decomposition.getLocalBoxMax();

  // Helper lambda to write single data values.
  auto printDataArray = [&](const auto &data, const std::string &type, const std::string &name) {
    vtsFile << "        <DataArray type=\"" << type << "\" Name=\"" << name << "\" format=\"ascii\">\n";
    vtsFile << "          " << data << "\n";
    vtsFile << "        </DataArray>\n";
  };

  vtsFile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n";
  vtsFile << "<VTKFile byte_order=\"LittleEndian\" type=\"StructuredGrid\" version=\"0.1\">\n";
  vtsFile << "  <StructuredGrid WholeExtent=\""
          << "0"
          << " " << dimensions[0] << " "
          << "0"
          << " " << dimensions[1] << " "
          << "0"
          << " " << dimensions[2] << "\">\n";
  vtsFile << "    <Piece Extent=\"" << (rankGridCoords[0]) << " " << (rankGridCoords[0] + 1) << " "
          << (rankGridCoords[1]) << " " << (rankGridCoords[1] + 1) << " " << (rankGridCoords[2]) << " "
          << (rankGridCoords[2] + 1) << "\">\n";
  vtsFile << "      <CellData>\n";
  printDataArray(autoPasConfig.cellSizeFactor, "Float32", "CellSizeFactor");
  printDataArray(static_cast<int>(autoPasConfig.container), "Int32", "Container");
  printDataArray(static_cast<int>(autoPasConfig.dataLayout), "Int32", "DataLayout");
  printDataArray(static_cast<int>(autoPasConfig.loadEstimator), "Int32", "LoadEstimator");
  printDataArray(static_cast<int>(autoPasConfig.traversal), "Int32", "Traversal");
  printDataArray(static_cast<int>(autoPasConfig.newton3), "Int32", "Newton3");
  printDataArray(rank, "Int32", "Rank");
  vtsFile << "      </CellData>\n";
  vtsFile << "      <Points>\n";
  vtsFile << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  vtsFile << "          " << localBoxMin[0] << " " << localBoxMin[1] << " " << localBoxMin[2] << "\n";
  vtsFile << "          " << localBoxMax[0] << " " << localBoxMin[1] << " " << localBoxMin[2] << "\n";
  vtsFile << "          " << localBoxMin[0] << " " << localBoxMax[1] << " " << localBoxMin[2] << "\n";
  vtsFile << "          " << localBoxMax[0] << " " << localBoxMax[1] << " " << localBoxMin[2] << "\n";
  vtsFile << "          " << localBoxMin[0] << " " << localBoxMin[1] << " " << localBoxMax[2] << "\n";
  vtsFile << "          " << localBoxMax[0] << " " << localBoxMin[1] << " " << localBoxMax[2] << "\n";
  vtsFile << "          " << localBoxMin[0] << " " << localBoxMax[1] << " " << localBoxMax[2] << "\n";
  vtsFile << "          " << localBoxMax[0] << " " << localBoxMax[1] << " " << localBoxMax[2] << "\n";
  vtsFile << "        </DataArray>\n";
  vtsFile << "      </Points>\n";
  vtsFile << "    </Piece>\n";
  vtsFile << "  </StructuredGrid>\n";
  vtsFile << "</VTKFile>\n";

  vtsFile.close();
}

std::string LADDS::RegularGridDecompositionLogger::filenameMetadata(size_t iteration) const {
  std::stringstream ss;
  ss << "Decomp_" << std::setfill('0') << std::setw(static_cast<int>(maxDigitsIterations)) << iteration << ".pvts";
  return ss.str();
}

std::string LADDS::RegularGridDecompositionLogger::filenamePayload(size_t iteration, int rank) const {
  if (rank == -1) {
    autopas::AutoPas_MPI_Comm_rank(decomposition.getCommunicator(), &rank);
  }
  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(decomposition.getCommunicator(), &numRanks);
  std::stringstream ss;
  ss << "Decomp_" << std::setfill('0') << std::setw(static_cast<int>(std::to_string(numRanks).size())) << rank << "_" 
     << std::setw(static_cast<int>(maxDigitsIterations)) << iteration << ".vts";
  return ss.str();
}
