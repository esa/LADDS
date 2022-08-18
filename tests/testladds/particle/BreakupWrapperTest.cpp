/**
 * @file BreakupWrapperTest.cpp
 * @author F. Gratl
 * @date 24.01.22
 */

#include "BreakupWrapperTest.h"

#include <gmock/gmock-matchers.h>

#include "ladds/particle/BreakupWrapper.h"

/**
 * Crash two particles into each other and observe that new particles have higher IDs.
 */
TEST_F(BreakupWrapperTest, testSimulationLoop) {
  auto [csvWriter, accumulator, integrator] = simulation.initIntegrator(*autopas, *configReader);

  // two particles 1000km above earth whose paths cross exactly at [R+1000, 0, 0]
  size_t highestIdBeforeCrash = 0;
  autopas->addParticle(LADDS::Particle({Physics::R_EARTH + 1000., -1., 0.},
                                       {0., 2., 0.},
                                       highestIdBeforeCrash,
                                       "A",
                                       LADDS::Particle::ActivityState::passive,
                                       1.,
                                       1.,
                                       LADDS::Particle::calculateBcInv(0., 1., 1., 2.2),
                                       std::numeric_limits<size_t>::max()));
  autopas->addParticle(LADDS::Particle({Physics::R_EARTH + 1000., 0., -1.},
                                       {0., 0., 2.},
                                       ++highestIdBeforeCrash,
                                       "B",
                                       LADDS::Particle::ActivityState::passive,
                                       1.,
                                       1.,
                                       LADDS::Particle::calculateBcInv(0., 1., 1., 2.2),
                                       std::numeric_limits<size_t>::max()));

  // dummy
  std::vector<LADDS::Constellation> constellations;
  // do one loop where we expect a breakup to happen
  simulation.simulationLoop(*autopas, *integrator, constellations, *configReader, *decomposition);

  // expect more particles than before
  EXPECT_GT(autopas->getNumberOfParticles(), 2);
  // expect every new particle to have a higher id than the particles before. Old particles are destroyed in the crash
  for (const auto &p : *autopas) {
    EXPECT_GT(p.getID(), highestIdBeforeCrash);
  }
}

TEST_F(BreakupWrapperTest, stressTest) {
  const size_t numCollisions = 100;
  const size_t numParticles = 2 * numCollisions;
  for (int i = 0; i < numParticles; ++i) {
    autopas->addParticle(LADDS::Particle({Physics::R_EARTH + 1000., 0., static_cast<double>(i)},
                                         {0., 2., 0.},
                                         i,
                                         "A",
                                         LADDS::Particle::ActivityState::passive,
                                         10.,
                                         1.,
                                         LADDS::Particle::calculateBcInv(0., 1., 1., 2.2),
                                         std::numeric_limits<size_t>::max()));
  }
  ASSERT_EQ(numParticles, autopas->getNumberOfParticles());
  std::vector<LADDS::Particle *> particlePointers;
  particlePointers.reserve(numParticles);
  autopas->forEach([&](auto &p) { particlePointers.push_back(&p); });

  LADDS::BreakupWrapper breakupWrapper(*configReader, *autopas);
  EXPECT_NO_THROW(for (int collisionId = 0; collisionId < numCollisions; ++collisionId) {
    breakupWrapper.simulateBreakup(
        {{particlePointers[collisionId], particlePointers[numCollisions / 2 + collisionId], 4.2, {0., 0., 0.}}});
  });
}
