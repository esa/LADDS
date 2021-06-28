/**
 * @file CollisionFunctor.h
 * @author F. Gratl
 * @date 28.06.21
 */

#include <autopas/AutoPasDecl.h>

#include <iostream>

#include "CollisionFunctor.h"
#include "Debris.h"
#include "Logger.h"

// Declare the main AutoPas class as extern template instantiation. It is instantiated in AutoPasClass.cpp.
extern template class autopas::AutoPas<Debris>;

int main() {
  Logger logger;

  // initialization of the simulation setup
  // TODO Read input
  constexpr size_t numDebris = 2;
  constexpr double cutoff = 2;

  // initialization of autopas
  autopas::AutoPas<Debris> autopas;
  autopas.setBoxMin({0., 0., 0.});
  autopas.setBoxMax({10., 10., 10.});
  autopas.setCutoff(cutoff);
  autopas.init();

  // initialization of the scenario
  for (size_t i = 0; i < numDebris; ++i) {
    autopas.addParticle(Debris{{static_cast<double>(i), 0, 0}, {0., 0., 0.}, i});
  }

  // just for fun: print particles
  for (const auto &d : autopas) {
    logger.log(Logger::Level::info, d.toString());
  }

  // pairwise interaction
  CollisionFunctor collisionFunctor(cutoff);
  autopas.iteratePairwise(&collisionFunctor);

  // print result
  logger.log(Logger::Level::info, "Close encounters: {}", collisionFunctor.getCollisions().size());

  return 0;
}
