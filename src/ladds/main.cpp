/**
 * @file CollisionFunctor.h
 * @author F. Gratl
 * @date 28.06.21
 */

#include <autopas/AutoPasDecl.h>
#include <breakupModel/output/VTKWriter.h>
#include <satellitePropagator/io/FileOutput.h>
#include <satellitePropagator/physics/AccelerationAccumulator.h>
#include <satellitePropagator/physics/Integrator.h>

#include <iostream>

#include "CollisionFunctor.h"
#include "DatasetReader.h"
#include "LoadConfig.h"
#include "Logger.h"
#include "Particle.h"
#include "SatelliteToParticleConverter.h"
#include "spdlog/fmt/ostr.h"

// Declare the main AutoPas class as extern template instantiation. It is instantiated in AutoPasClass.cpp.
extern template class autopas::AutoPas<Particle>;

int main(int argc, char **argv) {
  Logger logger;
  logger.get()->set_level(spdlog::level::debug);

  // Default config path
  auto cfgFilePath = LoadConfig::defaultCfgPath;

  // Read in config if given
  if (argc > 1) cfgFilePath = argv[1];
  const auto config = LoadConfig::loadConfig(cfgFilePath, logger);
  logger.log(Logger::Level::info, "Config loaded.");

  const size_t iterations = config["sim"]["iterations"].as<size_t>();
  const size_t vtkWriteFrequency = config["io"]["vtkWriteFrequency"].as<size_t>();

  using AutoPas_t = autopas::AutoPas<Particle>;

  // initialization of autopas
  AutoPas_t autopas;
  const double maxAltitude = config["sim"]["maxAltitude"].as<double>();
  const double cutoff = config["autopas"]["cutoff"].as<double>();
  const double desiredCellsPerDimension = config["autopas"]["desiredCellsPerDimension"].as<double>();

  autopas.setBoxMin({-maxAltitude, -maxAltitude, -maxAltitude});
  autopas.setBoxMax({maxAltitude, maxAltitude, maxAltitude});
  autopas.setCutoff(cutoff);
  // Set the size (relative to cutoff) of the cells so that roughly the desired number of cells per dimension is reached
  autopas.setCellSizeFactor((maxAltitude * 2.) / (cutoff * desiredCellsPerDimension));
  autopas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autopas.setAllowedContainers({autopas::ContainerOption::verletListsCells});
  autopas.init();

  // initialization of the integrator
  std::array<bool, 8> selectedPropagatorComponents{
      config["sim"]["prop"]["useKEPComponent"].as<bool>(), config["sim"]["prop"]["useJ2Component"].as<bool>(),
      config["sim"]["prop"]["useC22Component"].as<bool>(), config["sim"]["prop"]["useS22Component"].as<bool>(),
      config["sim"]["prop"]["useSOLComponent"].as<bool>(), config["sim"]["prop"]["useLUNComponent"].as<bool>(),
      config["sim"]["prop"]["useSRPComponent"].as<bool>(), config["sim"]["prop"]["useDRAGComponent"].as<bool>()};

  auto fo = std::make_shared<FileOutput<AutoPas_t, Particle>>(autopas, config["io"]["output_file"].as<std::string>(),
                                                              OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t, Particle>>(
      selectedPropagatorComponents, autopas, 0.0, *fo);
  auto integrator = std::make_shared<Integrator<AutoPas_t, Particle>>(autopas, *accumulator, 1e-1);

  // Read in scenario
  auto actualSatellites =
      DatasetReader::readDataset(std::string(DATADIR) + config["io"]["posFileName"].as<std::string>(),
                                 std::string(DATADIR) + config["io"]["velFileName"].as<std::string>());
  logger.log(Logger::Level::debug, "Parsed {} satellites", actualSatellites.size());

  double minAltitudeFound{std::numeric_limits<double>::max()};
  double maxAltitudeFound{0.};
  // Convert satellites to particles
  for (const auto &particle : actualSatellites) {
    auto pos = particle.getPosition();
    double altitude = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
    minAltitudeFound = std::max(minAltitudeFound, altitude);
    maxAltitudeFound = std::max(maxAltitudeFound, altitude);
    if (altitude < maxAltitude) {
      autopas.addParticle(particle);
    }
  }
  logger.log(Logger::Level::info, "Min altitude is {}", minAltitudeFound);
  logger.log(Logger::Level::info, "Max altitude is {}", maxAltitudeFound);

  logger.log(Logger::Level::info, "Number of particles: {}", autopas.getNumberOfParticles());

  // main-loop skeleton
  for (size_t i = 0; i < iterations; ++i) {
    // update positions
    integrator->integrate(false);

    // TODO MPI: handle particle exchange between ranks
    // (potentially) update the internal data structure and collect particles which are leaving the container.
    const auto [escapedParticles, containerUpdated] = autopas.updateContainer();

    if (not escapedParticles.empty()) {
      logger.log(Logger::Level::err, "Particles are escaping! \n{}", escapedParticles);
    }
    // TODO Check for particles that burn up
    // pairwise interaction
    CollisionFunctor collisionFunctor(cutoff);
    autopas.iteratePairwise(&collisionFunctor);
    auto collisions = collisionFunctor.getCollisions();
    logger.log(Logger::Level::info, "Close encounters: {}", collisions.size());
    for (const auto &[p1, p2] : collisions) {
      logger.log(Logger::Level::debug, "{} | {}", p1->getID(), p2->getID());
    }

    // TODO insert breakup model here

    // Visualization:
    if (i % vtkWriteFrequency == 0) {
      VTKWriter vtkWriter("output_" + std::to_string(i) + ".vtu");
      std::vector<Satellite> allParticles;
      allParticles.reserve(autopas.getNumberOfParticles());
      for (const auto &p : autopas) {
        allParticles.push_back(SatelliteToParticleConverter::convertParticleToSatellite(p));
      }
      vtkWriter.printResult(allParticles);
    }
  }

  return 0;
}
