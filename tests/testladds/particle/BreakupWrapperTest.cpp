/**
 * @file BreakupWrapperTest.cpp
 * @author F. Gratl
 * @date 24.01.22
 */

#include "BreakupWrapperTest.h"

#include <gmock/gmock-matchers.h>

/**
 * Crash two particles into each other and observe that new particles have higher IDs.
 */
TEST_F(BreakupWrapperTest, testSimulationLoop) {
  auto [csvWriter, accumulator, integrator] = simulation.initIntegrator(*autopas, *configReader);

  // two particles 1000km above earth whose paths cross exactly at [R+1000, 0, 0]
  size_t highestIdBeforeCrash = 4;
  autopas->addParticle(LADDS::Particle({Physics::R_EARTH + 1000., -1., 0.},
                                       {0., 2., 0.},
                                       1,
                                       "A",
                                       LADDS::Particle::ActivityState::passive,
                                       1.,
                                       1.,
                                       LADDS::Particle::calculateBcInv(0., 1., 1., 2.2)));
  autopas->addParticle(LADDS::Particle({Physics::R_EARTH + 1000., 0., -1.},
                                       {0., 0., 2.},
                                       highestIdBeforeCrash,
                                       "B",
                                       LADDS::Particle::ActivityState::passive,
                                       1.,
                                       1.,
                                       LADDS::Particle::calculateBcInv(0., 1., 1., 2.2)));

  // dummy
  std::vector<LADDS::Constellation> constellations;
  // do one loop where we expect a breakup to happen
  simulation.simulationLoop(*autopas, *integrator, constellations, *configReader);

  // expect more particles than before
  EXPECT_GT(autopas->getNumberOfParticles(), 2);
  // expect every new particle to have a higher id than the particles before. Old particles are destroyed in the crash
  for (const auto &p : *autopas) {
    EXPECT_GT(p.getID(), highestIdBeforeCrash);
  }
}
