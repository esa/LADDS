/**
 * @file SatelliteLoader.cpp
 * @author F. Gratl
 * @date 06.12.21
 */

#include "SatelliteLoader.h"

#include "ConfigReader.h"
#include "DatasetReader.h"
#include "Logger.h"
#include "ladds/distributedMemParallelization/AltitudeBasedDecomposition.h"
#include "ladds/io/hdf5/HDF5Reader.h"

namespace LADDS {

void SatelliteLoader::addSatellitesToAutoPas(AutoPas_t &autopas,
                                             std::vector<Particle> &satellites,
                                             DomainDecomposition &decomp,
                                             ConfigReader &config) {
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);
  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numRanks);

  // We rebalance our domain decomposition based on particle positions. (only makes sense with enough particles)
  if (satellites.size() >= (unsigned)numRanks) {
    decomp.rebalanceDecomposition(satellites, autopas);
  }

  // load particle vector into autopas while checking that they are within the desired altitude
  const auto maxAltitude = config.get<double>("sim/maxAltitude");
  const auto maxAltitudeSquared = maxAltitude * maxAltitude;
  double minAltitudeFound{std::numeric_limits<double>::max()};
  double maxAltitudeFound{0.};
  const auto &boxMin = autopas.getBoxMin();
  const auto &boxMax = autopas.getBoxMax();

  SPDLOG_LOGGER_INFO(config.getLogger().get(), "BoxMin: {}, BoxMax: {}", boxMin, boxMax);

  for (auto &particle : satellites) {
    const auto &pos = particle.getPosition();
    const auto rankForParticle = decomp.getRank(pos);
    if (rankForParticle != rank) {
      continue;
    }

    const double altitudeSquared = autopas::utils::ArrayMath::dot(pos, pos);
    minAltitudeFound = std::min(minAltitudeFound, altitudeSquared);
    maxAltitudeFound = std::max(maxAltitudeFound, altitudeSquared);

    if (altitudeSquared < maxAltitudeSquared) {
      // set the id to the next free particle id of this rank
      particle.setID(config.newParticleID());
      autopas.addParticle(particle);
    } else {
      SPDLOG_LOGGER_WARN(config.getLogger().get(),
                         "Particle NOT added because its altitudeSquared was too high!\n"
                         "Max allowed: {}\n"
                         "Actual: {}\n"
                         "{})",
                         maxAltitude,
                         sqrt(altitudeSquared),
                         particle.toString());
    }
  }

  // collect and reduce global values
  double minAltitudeFoundGlobal{};
  autopas::AutoPas_MPI_Reduce(
      &minAltitudeFound, &minAltitudeFoundGlobal, 1, AUTOPAS_MPI_DOUBLE, AUTOPAS_MPI_MIN, 0, decomp.getCommunicator());
  double maxAltitudeFoundGlobal{};
  autopas::AutoPas_MPI_Reduce(
      &maxAltitudeFound, &maxAltitudeFoundGlobal, 1, AUTOPAS_MPI_DOUBLE, AUTOPAS_MPI_MAX, 0, decomp.getCommunicator());
  minAltitudeFoundGlobal = sqrt(minAltitudeFoundGlobal);
  maxAltitudeFoundGlobal = sqrt(maxAltitudeFoundGlobal);
  SPDLOG_LOGGER_INFO(config.getLogger().get(), "Min altitude is {}", minAltitudeFoundGlobal);
  SPDLOG_LOGGER_INFO(config.getLogger().get(), "Max altitude is {}", maxAltitudeFoundGlobal);
  unsigned long numParticles = autopas.getNumberOfParticles();
  unsigned long numParticlesGlobal{};
  autopas::AutoPas_MPI_Reduce(
      &numParticles, &numParticlesGlobal, 1, AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_SUM, 0, decomp.getCommunicator());
  SPDLOG_LOGGER_INFO(config.getLogger().get(), "Number of particles: {}", numParticlesGlobal);
}

void SatelliteLoader::loadSatellites(AutoPas_t &autopas, ConfigReader &config, DomainDecomposition &decomp) {
  std::vector<Particle> satellites;

  // if the value is not set this is the only place that should define the default value as it is the first to access it
  const auto coefficientOfDrag = config.get<double>("sim/prop/coefficientOfDrag", 2.2);

  const auto csvFileName = config.get<std::string>("io/csv/fileName", "", true);
  // avoid creating io/hdf5 in the config if it doesn't exist
  const auto checkpointPathCfg =
      config.defines("io/hdf5", true) ? config.get<std::string>("io/hdf5/checkpoint/file", "", true) : "";
  // sanity checks
  if (csvFileName.empty() and checkpointPathCfg.empty()) {
    throw std::runtime_error("No valid input option found! Exiting...");
  }
  if (not csvFileName.empty() and not checkpointPathCfg.empty()) {
    throw std::runtime_error("Ambiguous input option found! Both CSV and HDF5 checkpoint found. Exiting...");
  }

  // load CSV ...
  if (not csvFileName.empty()) {
    // TODO: Remove DATADIR functionality
    const auto csvFilePath = std::string(DATADIR) + csvFileName;

    SPDLOG_LOGGER_INFO(config.getLogger().get(), "Loading scenario from CSV: {}", csvFilePath);

    // Read in scenario
    satellites = DatasetReader::readDataset(csvFilePath, coefficientOfDrag);
    SPDLOG_LOGGER_DEBUG(config.getLogger().get(), "Parsed {} satellites", satellites.size());
  }
  // ... or load checkpoint ...
  if (not checkpointPathCfg.empty()) {
    // TODO: Remove DATADIR functionality
    const auto checkpointPath = std::string(DATADIR) + checkpointPathCfg;
    SPDLOG_LOGGER_INFO(config.getLogger().get(), "Loading scenario from HDF5 checkpoint\nFile: {}", checkpointPath);

    HDF5Reader hdfReader(checkpointPath);
    // either load the given iteration or fall back to the last iteration stored in the file
    auto iteration = config.get<size_t>("io/hdf5/checkpoint/iteration", hdfReader.readLastIterationNr());
    satellites = hdfReader.readParticles(iteration, coefficientOfDrag);
  }

  addSatellitesToAutoPas(autopas, satellites, decomp, config);
}

std::vector<Constellation> SatelliteLoader::loadConstellations(ConfigReader &config, const Logger &logger) {
  std::vector<Constellation> constellations;
  auto constellationList = config.get<std::string>("io/constellationList", "", true);
  if (not constellationList.empty()) {
    auto constellationDataStr = config.get<std::string>("io/constellationList");
    // count constellation by counting ';'
    int nConstellations = std::count(constellationDataStr.begin(), constellationDataStr.end(), ';') + 1;

    // parse constellation info
    constellations.reserve(nConstellations);
    for (int i = 0; i < nConstellations; ++i) {
      unsigned long offset =
          (i == nConstellations - 1) ? constellationDataStr.size() : constellationDataStr.find(';', 0);
      std::string constellationDir = constellationDataStr.substr(0, offset);

      ConfigReader constellationConfig =
          ConfigReader(std::string(DATADIR) + constellationDir + "/shells_" + constellationDir + ".yaml", logger);

      constellations.emplace_back(Constellation(constellationConfig, config));
      if (i != nConstellations - 1) {
        constellationDataStr.erase(0, offset + 1);
      }
    }

    size_t constellationTotalNumSatellites = 0;
    for (const auto &constellation : constellations) {
      constellationTotalNumSatellites += constellation.getConstellationSize();
    }

    SPDLOG_LOGGER_INFO(logger.get(),
                       "{} more particles will be added from {} constellations",
                       constellationTotalNumSatellites,
                       nConstellations);

    for (const auto &c : constellations) {
      SPDLOG_LOGGER_INFO(
          logger.get(),
          "{}: insertion starts at iteration: {}, is fully deployed within {} iterations, inserts {} satellites",
          c.getConstellationName(),
          c.getStartTime(),
          c.getDuration(),
          c.getConstellationSize());
    }
  }
  return constellations;
}

}  // namespace LADDS