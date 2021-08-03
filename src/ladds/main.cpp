/**
 * @file CollisionFunctor.h
 * @author F. Gratl
 * @date 28.06.21
 */

#include <autopas/AutoPasDecl.h>
#include <breakupModel/input/TLESatcatDataReader.h>
#include <breakupModel/output/VTKWriter.h>
#include <satellitePropagator/io/FileOutput.h>
#include <satellitePropagator/physics/AccelerationAccumulator.h>
#include <satellitePropagator/physics/Integrator.h>

#include <iostream>

#include "CollisionFunctor.h"
#include "Logger.h"
#include "Particle.h"
#include "SatelliteToParticleConverter.h"
#include "spdlog/fmt/ostr.h"

// Declare the main AutoPas class as extern template instantiation. It is instantiated in AutoPasClass.cpp.
extern template class autopas::AutoPas<Particle>;

int main(int argc, char **argv) {
  Logger logger;
  logger.get()->set_level(spdlog::level::debug);

  // initialization of the simulation setup
  // TODO Read input
  std::string tleFilePath;
  if(argc > 1) {
    tleFilePath = argv[1];
  }
  if(tleFilePath.empty()){
    logger.log(Logger::Level::critical, "No TLE file given!");
  }
  constexpr double cutoff = 0.02;
  const size_t iterations = 3;
  const double max_altitude = 85000.;
  const double desiredCellsPerDimension = 50;
  const size_t vtkWriteFrequency = 10;

  using AutoPas_t = autopas::AutoPas<Particle>;

  // initialization of autopas
  AutoPas_t autopas;
  autopas.setBoxMin({-max_altitude, -max_altitude, -max_altitude});
  autopas.setBoxMax({max_altitude, max_altitude, max_altitude});
  autopas.setCutoff(cutoff);
  autopas.setCellSizeFactor((max_altitude * 2.) / (cutoff * desiredCellsPerDimension));
  autopas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autopas.setAllowedContainers({autopas::ContainerOption::verletListsCells});
  autopas.init();

  // initialization of the integrator
  std::array<bool, 8> selectedPropagatorComponents{true, false, false, false, false, false, false, false};
  auto fo = std::make_shared<FileOutput<AutoPas_t, Particle>>(autopas, "output.csv", OutputFile::CSV,
                                                              selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t, Particle>>(
      selectedPropagatorComponents, autopas, 0.0, *fo);
  auto integrator = std::make_shared<Integrator<AutoPas_t, Particle>>(autopas, *accumulator, 1e-1);

  // Read in scenario
  TLESatcatDataReader tleSatcatDataReader{std::string(DATADIR) + "satcat_breakupModel.csv", tleFilePath};
  auto actualSatellites = tleSatcatDataReader.getSatelliteCollection();
  logger.log(Logger::Level::debug, "Parsed {} satellites", actualSatellites.size());

  double minAltitudeFound{std::numeric_limits<double>::max()};
  double maxAltitudeFound{0.};
  // Convert satellites to particles
  for (const auto &satellite : actualSatellites) {
    auto particle = SatelliteToParticleConverter::convertSatelliteToParticle(satellite);
    auto pos = particle.getPosition();
    double altitude = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
    minAltitudeFound = std::max(minAltitudeFound, altitude);
    maxAltitudeFound = std::max(maxAltitudeFound, altitude);
    if (altitude < max_altitude) {
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
    const auto [escapedParticles, containerUpdated] = autopas.updateContainer();

    if (not escapedParticles.empty()) {
      logger.log(Logger::Level::err, "Particles are escaping! \n{}", escapedParticles);
    }
//TODO Check for particles that burn up
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
