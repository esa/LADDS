/**
 * @file CollisionFunctor.h
 * @author F. Gratl
 * @date 28.06.21
 */

#include <autopas/utils/ArrayUtils.h>
#include <satellitePropagator/physics/Integrator.h>
#include <satellitePropagator/physics/AccelerationAccumulator.h>
#include <satellitePropagator/io/FileOutput.h>

#include <breakupModel/input/TLESatcatDataReader.h>


#include <autopas/AutoPasDecl.h>

#include <iostream>

#include "CollisionFunctor.h"
#include "Particle.h"
#include "Logger.h"
#include "SatelliteToParticleConverter.h"

// Declare the main AutoPas class as extern template instantiation. It is instantiated in AutoPasClass.cpp.
extern template class autopas::AutoPas<Particle>;

int main() {
  Logger logger;

  // initialization of the simulation setup
  // TODO Read input
  constexpr size_t numDebris = 2;
  constexpr double cutoff = 0.02;
  const size_t iterations = 1000;
  const double max_altitude = 8500.;
  double min_altitude = 1.e9;

  using AutoPas_t = autopas::AutoPas<Particle>;

  // initialization of autopas
  AutoPas_t autopas;
  autopas.setBoxMin({-max_altitude, -max_altitude, -max_altitude});
  autopas.setBoxMax({max_altitude, max_altitude, max_altitude});
  autopas.setCutoff(cutoff);
  autopas.setCellSizeFactor(20000);
  autopas.init();

  // initialization of the integrator
  std::array<bool,8> selectedPropagatorComponents{true, false, false, false, false, false, false, false};
  auto fo = std::make_shared<FileOutput<AutoPas_t,Particle>>(autopas, "output.csv",OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t,Particle>>(selectedPropagatorComponents,autopas,0.0,*fo);
  auto integrator = std::make_shared<Integrator<AutoPas_t,Particle>>(autopas,*accumulator,1e-4);


  // Read in scenario
  TLESatcatDataReader tleSatcatDataReader{"/home/pablo/Code/LADDS/build/_deps/breakupmodelfetch-src/satcat.csv", "/home/pablo/Code/LADDS/data/tle.txt"};
  auto actualSatellites = tleSatcatDataReader.getSatelliteCollection();

  // Convert satellites to particles
  for (const auto &satellite : actualSatellites) {
    auto particle = SatelliteToParticleConverter::convertSatelliteToParticle(satellite);
    auto pos = particle.getPosition();
    double altitude = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    if (altitude < min_altitude) {
      min_altitude = altitude;
    }
    if (altitude < max_altitude) {
      autopas.addParticle(particle);
    }
  }
  logger.log(Logger::Level::info, "Min_altitude was {}", min_altitude);

  // just for fun: print particles
  // for (const auto &d : autopas) {
  //   logger.log(Logger::Level::info, d.toString());
  // }
  logger.log(Logger::Level::info, "Number of particles: {}", autopas.getNumberOfParticles());

  // main-loop skeleton
  for (size_t i = 0; i < iterations; ++i) {
    // update positions
    integrator->integrate(false);

    std::cout << autopas::utils::ArrayUtils::to_string(autopas.begin()->getPosition()) << std::endl;

    // TODO MPI: handle particle exchange between ranks

    // pairwise interaction
    CollisionFunctor collisionFunctor(cutoff);
    autopas.iteratePairwise(&collisionFunctor);
    logger.log(Logger::Level::info, "Close encounters: {}", collisionFunctor.getCollisions().size());

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

  // just for fun: print particles
  // for (const auto &d : autopas) {
  //   logger.log(Logger::Level::info, d.toString());
  // }

  return 0;
}
