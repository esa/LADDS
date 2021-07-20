/**
 * @file CollisionFunctor.h
 * @author F. Gratl
 * @date 28.06.21
 */

#include <autopas/AutoPasDecl.h>

#include <iostream>

#include "CollisionFunctor.h"
#include "Particle.h"
#include "Logger.h"

#include <satellitePropagator/physics/Integrator.h>
#include <satellitePropagator/physics/AccelerationAccumulator.h>
#include <satellitePropagator/io/FileOutput.h>

// Declare the main AutoPas class as extern template instantiation. It is instantiated in AutoPasClass.cpp.
extern template class autopas::AutoPas<Particle>;

int main() {
  Logger logger;

  // initialization of the simulation setup
  // TODO Read input
  constexpr size_t numDebris = 2;
  constexpr double cutoff = 2;
  const size_t iterations = 1;

  using AutoPas_t = autopas::AutoPas<Particle>;

  // initialization of autopas
  AutoPas_t autopas;
  autopas.setBoxMin({0., 0., 0.});
  autopas.setBoxMax({10., 10., 10.});
  autopas.setCutoff(cutoff);
  autopas.init();

  // initialization of the integrator
  std::array<bool,8> selectedPropagatorComponents{true, false, false, false, false, false, false, false};
  auto fo = std::make_shared<FileOutput<AutoPas_t,Particle>>(autopas, "output.csv",OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t,Particle>>(selectedPropagatorComponents,autopas,0.0,*fo);
  auto integrator = std::make_shared<Integrator<AutoPas_t,Particle>>(autopas,*accumulator,1e-8);


  // initialization of the scenario
  for (size_t i = 0; i < numDebris; ++i) {
    autopas.addParticle(Particle{{static_cast<double>(i), 0, 0}, {0., 0., 0.}, i});
  }

  // just for fun: print particles
  for (const auto &d : autopas) {
    logger.log(Logger::Level::info, d.toString());
  }

  // main-loop skeleton
  for (size_t i = 0; i < iterations; ++i) {
    // update positions
    // TODO insert propagator code here
    // TODO MPI: handle particle exchange between ranks

    // pairwise interaction
    CollisionFunctor collisionFunctor(cutoff);
    autopas.iteratePairwise(&collisionFunctor);
    logger.log(Logger::Level::info, "Close encounters: {}", collisionFunctor.getCollisions().size());

    // TODO insert breakup model here
  }

  return 0;
}
