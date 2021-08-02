/**
 * @file CollisionFunctor.h
 * @author F. Gratl
 * @date 28.06.21
 */

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
  constexpr double cutoff = 500.;
  const size_t iterations = 1;

  using AutoPas_t = autopas::AutoPas<Particle>;

  // initialization of autopas
  AutoPas_t autopas;
  autopas.setBoxMin({0., 0., 0.});
  autopas.setBoxMax({10000., 10000., 10000.});
  autopas.setCutoff(cutoff);
  autopas.init();

  // initialization of the integrator
  std::array<bool,8> selectedPropagatorComponents{true, false, false, false, false, false, false, false};
  auto fo = std::make_shared<FileOutput<AutoPas_t,Particle>>(autopas, "output.csv",OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t,Particle>>(selectedPropagatorComponents,autopas,0.0,*fo);
  auto integrator = std::make_shared<Integrator<AutoPas_t,Particle>>(autopas,*accumulator,1e-8);


  // Read in scenario
  TLESatcatDataReader tleSatcatDataReader{"/home/pablo/Code/LADDS/build/_deps/breakupmodelfetch-src/satcat.csv", "/home/pablo/Code/LADDS/data/tle.txt"};
  auto actualSatellites = tleSatcatDataReader.getSatelliteCollection();

  // Convert satellites to particles
  for (const auto &satellite : actualSatellites) {
    autopas.addParticle(SatelliteToParticleConverter::convertSatelliteToParticle(satellite));
  }

  // just for fun: print particles
  for (const auto &d : autopas) {
    logger.log(Logger::Level::info, d.toString());
  }

  // main-loop skeleton
  for (size_t i = 0; i < iterations; ++i) {
    // update positions
    integrator->integrate(false);

    // TODO MPI: handle particle exchange between ranks

    // pairwise interaction
    CollisionFunctor collisionFunctor(cutoff);
    autopas.iteratePairwise(&collisionFunctor);
    logger.log(Logger::Level::info, "Close encounters: {}", collisionFunctor.getCollisions().size());

    // TODO insert breakup model here
  }

  // just for fun: print particles
  for (const auto &d : autopas) {
    logger.log(Logger::Level::info, d.toString());
  }

  return 0;
}
