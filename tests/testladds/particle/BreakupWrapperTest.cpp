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

  // two particles whose paths cross exactly
  size_t highestIdBeforeCrash = 4;
  autopas->addParticle(Particle({Physics::R_EARTH + 1000., -1., 0.}, {0., 2., 0.}, 1));
  autopas->addParticle(Particle({Physics::R_EARTH + 1000., 0., -1.}, {0., 0., 2.}, highestIdBeforeCrash));

  // dummy
  std::vector<Constellation> constellations;
  // do one loop where we expect a breakup to happen
  simulation.simulationLoop(*autopas, *integrator, constellations, *configReader);

  // expect more particles than before
  EXPECT_GT(autopas->getNumberOfParticles(), 2);
  // expect every new particle to have a higher id than the particles before. Old particles are destroyed in the crash
  for (const auto &p : *autopas) {
    EXPECT_GT(p.getID(), highestIdBeforeCrash);
  }
}