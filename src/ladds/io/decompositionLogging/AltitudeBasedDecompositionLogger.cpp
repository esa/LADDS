/**
 * @file AltitudeBasedDecompositionLogger.cpp
 * @author F. Gratl
 * @date 24.10.22
 */

#include "AltitudeBasedDecompositionLogger.h"

LADDS::AltitudeBasedDecompositionLogger::AltitudeBasedDecompositionLogger(
    LADDS::ConfigReader &config, const LADDS::AltitudeBasedDecomposition &decomposition)
    : DecompositionLoggerParametrized<AltitudeBasedDecomposition>(
          "AltitudeBasedDecompositionLogger", config, "vtp", decomposition) {}

void LADDS::AltitudeBasedDecompositionLogger::writeMetafile(size_t iteration) const {
  const auto filename = filenameMetadata(iteration);

  std::ofstream pvtsFile;
  pvtsFile.open(filename, std::ios::out);

  if (not pvtsFile.is_open()) {
    throw std::runtime_error("AltitudeBasedDecompositionLogger::writeMetafile(): Failed to open file \"" + filename +
                             "\"");
  }

  const std::array<double, 3> globalBoxMin = decomposition.getGlobalBoxMin();
  const std::array<double, 3> globalBoxMax = decomposition.getGlobalBoxMax();
  pvtsFile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n";
  pvtsFile << "<VTKFile byte_order=\"LittleEndian\" type=\"PPolyData\" version=\"0.1\">\n";
  pvtsFile << "  <PPolyData GhostLevel=\"1\">\n";
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
    pvtsFile << "    <Piece Source=\"" << filenamePiece << "\"/>\n";
  }

  pvtsFile << "  </PPolyData>\n";
  pvtsFile << "</VTKFile>\n";

  pvtsFile.close();
}

void LADDS::AltitudeBasedDecompositionLogger::writePayload(size_t iteration,
                                                           const autopas::Configuration &autoPasConfig) const {
  const auto timestepFileName = filenamePayload(iteration);

  std::ofstream vtsFile;
  vtsFile.open(timestepFileName, std::ios::out);

  //  const auto &[dimensions, periods, rankGridCoords] = decomposition.getGridInfo();
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(decomposition.getCommunicator(), &rank);
  const std::array<double, 3> localBoxMin = decomposition.getLocalBoxMin();
  const std::array<double, 3> localBoxMax = decomposition.getLocalBoxMax();
  //  int numRanks{};
  //  autopas::AutoPas_MPI_Comm_size(decomposition.getCommunicator(), &numRanks);

  // Helper lambda to write single data values.
  auto printDataArray = [&](const auto &data, const std::string &type, const std::string &name) {
    vtsFile << "        <DataArray type=\"" << type << "\" Name=\"" << name << "\" format=\"ascii\">\n";
    vtsFile << "          " << data << " " << data << " " << data << " " << data << " " << data << " " << data << "\n";
    vtsFile << "        </DataArray>\n";
  };

  vtsFile << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n";
  vtsFile << "<VTKFile byte_order=\"LittleEndian\" type=\"PolyData\" version=\"0.1\">\n";
  vtsFile << "  <PolyData>\n";
  vtsFile << "    <Piece NumberOfPoints=\"" << 8
          << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\"\n"
             "NumberOfStrips=\"0\" NumberOfPolys=\"6\">\n";
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
  vtsFile << "      <Polys>\n";
  // which points are connected to one polygon, (read: 0->1, 0->1->3, ...)
  vtsFile << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  vtsFile << "          0 1 3 2 4 5 7 6 1 5 7 3 0 4 6 2 2 3 7 6 0 1 5 4\n";
  vtsFile << "        </DataArray>\n";
  // which entries from "connectivity form a face (read 1-4, 5-8, ...)
  vtsFile << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  vtsFile << "          4 8 12 16 20 24\n";
  vtsFile << "        </DataArray>\n";
  vtsFile << "      </Polys>\n";
  vtsFile << "    </Piece>\n";
  vtsFile << "  </PolyData>\n";
  vtsFile << "</VTKFile>\n";

  vtsFile.close();
}
