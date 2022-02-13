/**
 * @file SimulationTest.cpp
 * @author F. Gratl
 * @date 03.12.21
 */

#include "SimulationTest.h"

#include <autopas/AutoPasDecl.h>
#include <satellitePropagator/io/FileOutput.h>
#include <satellitePropagator/physics/AccelerationAccumulator.h>
#include <satellitePropagator/physics/Integrator.h>

#include "ladds/particle/Particle.h"

/**
 * Tests whether particles are correctly inserted into the simulation, when a particle
 * from the simulation and from the insertion overlap (or are very very close to each other)
 */
TEST_F(SimulationTest, testInsertionOverlap) {
  auto [csvWriter, accumulator, integrator] = simulation.initIntegrator(*autopas, *configReader);

  /*
   * position of p2 is the position of p1 after 1 second
   * p3,p4 are copies of p1,p2
   * -> the following formations are expected after iteration i:
   */

  /*  i |   sim   |  insertion
   *  0   -oo---    -oo----      expected simSize = 2
   *  1   --oo--    -oo----
   *  1'  -ooo--    --o----      expected simSize = 3
   *  2   --ooo-    --o----      expected simSize = 3
   *  3   ---ooo    --o----
   *  3'  --oooo    -------      expected simSize = 4
   *  4   ---oooo   -------      expected simSize = 4
   */
  const std::vector<size_t> expectedDelayedParticles{2, 1, 1, 0, 0, 0};

  const std::vector<Particle> initialSatellites{Particle({6871., 0., 0.},
                                                         {0, 4.8958309146899, 5.83462408131549},
                                                         1,
                                                         "initial 1",
                                                         Particle::ActivityState::passive,
                                                         1.,
                                                         1.,
                                                         Particle::calculateBcInv(0., 1., 1., 2.2)),
                                                Particle({6870.99577848984, 4.89582991243564, 5.83462288687537},
                                                         {-0.00844301944979238, 4.89582790671436, 5.83462049654983},
                                                         2,
                                                         "initial 2",
                                                         Particle::ActivityState::passive,
                                                         1.,
                                                         1.,
                                                         Particle::calculateBcInv(0., 1., 1., 2.2))};
  for (const auto &s : initialSatellites) {
    autopas->addParticle(s);
  }
  const std::vector<Particle> newSatellites{Particle({6871., 0., 0.},
                                                     {0, 4.8958309146899, 5.83462408131549},
                                                     1,
                                                     "new 1",
                                                     Particle::ActivityState::passive,
                                                     0,
                                                     0,
                                                     Particle::calculateBcInv(0., 1., 1., 2.2)),
                                            Particle({6870.99577848984, 4.89582991243564, 5.83462288687537},
                                                     {-0.00844301944979238, 4.89582790671436, 5.83462049654983},
                                                     4,
                                                     "new 2",
                                                     Particle::ActivityState::passive,
                                                     1.,
                                                     1.,
                                                     Particle::calculateBcInv(0., 1., 1., 2.2))};

  auto constellationCutoff = config["io"]["constellationCutoff"].as<double>();

  ASSERT_EQ(autopas->getNumberOfParticles(), 2) << "Container initialized wrong!";
  std::vector<Particle> delayedSatellites = newSatellites;
  // simulate some iterations and check delays against expectations from above
  for (size_t iteration = 0; iteration < expectedDelayedParticles.size(); ++iteration) {
    auto escapedParticles = autopas->updateContainer();
    ASSERT_EQ(escapedParticles.size(), 0) << "In this test nothing should escape.";
    delayedSatellites = simulation.checkedInsert(*autopas, delayedSatellites, constellationCutoff);
    EXPECT_EQ(delayedSatellites.size(), expectedDelayedParticles[iteration]);
    EXPECT_EQ(autopas->getNumberOfParticles(),
              initialSatellites.size() + newSatellites.size() - delayedSatellites.size());

    integrator->integrate();
  }
}

/**
 * Initialize two particles close to the burn up radius, one going towards earth one away from it.
 * Do one simulation iteration and check that exactly the right one was removed.
 */
TEST_F(SimulationTest, testBurnUp) {
  configReader->setValue("sim/breakup/enabled", true);
  auto [csvWriter, accumulator, integrator] = simulation.initIntegrator(*autopas, *configReader);
  const auto &minAltitude = Physics::R_EARTH + configReader->get<double>("sim/minAltitude");
  // initialize a particle 1km above burn up radius with a trajectory towards earth
  autopas->addParticle(Particle({minAltitude + 1., 0., 0.},
                                {-10., 0., 0.},
                                0,
                                "1km above ground towards earth",
                                Particle::ActivityState::passive,
                                1.,
                                1.,
                                Particle::calculateBcInv(0., 1., 1., 2.2)));
  // initialize a particle 1km above burn up radius with a trajectory away from earth
  // different position to avoid any interferences
  autopas->addParticle(Particle({0., minAltitude + 1., 0.},
                                {0., 10., 0.},
                                1,
                                "1km above ground away from earth",
                                Particle::ActivityState::passive,
                                1.,
                                1.,
                                Particle::calculateBcInv(0., 1., 1., 2.2)));
  ASSERT_EQ(autopas->getNumberOfParticles(), 2) << "Initial particles not found!";

  // dummy because interface requires it
  std::vector<Constellation> constellations{};
  // do one iteration where the burnup is expected to happen
  simulation.simulationLoop(*autopas, *integrator, constellations, *configReader);

  EXPECT_EQ(autopas->getNumberOfParticles(), 1) << "Exactly one particle should have burnt up!";
  EXPECT_EQ(autopas->begin()->getID(), 1) << "Remaining particle has unexpected Id!";
}

/**
 * Tests whether particles are only inserted, when they are outside a critical range
 * and delayed when within that critical range.
 *
 * For checked positions see generateParameters()
 */
TEST_P(SimulationTest, testCheckedInsert) {
  const auto &[posTestParticle, positionIsSafe] = GetParam();

  autopas->addParticle(Particle(testCheckedInsertParticlePos,
                                zeroVec,
                                0,
                                "existing",
                                Particle::ActivityState::passive,
                                1.,
                                1.,
                                Particle::calculateBcInv(0., 1., 1., 2.2)));

  // particle that will be inserted
  Particle p1{posTestParticle,
              zeroVec,
              1,
              "tester",
              Particle::ActivityState::passive,
              1.,
              1.,
              Particle::calculateBcInv(0., 1., 1., 2.2)};

  const auto escapedParticles = autopas->updateContainer();
  ASSERT_TRUE(escapedParticles.empty()) << "Test setup faulty!";

  std::vector<Particle> delayedParticles = simulation.checkedInsert(*autopas, {p1}, constellationCutoff);
  // if the position is safe the particle is inserted and there are now two. Otherwise, we remain with one.
  EXPECT_EQ(autopas->getNumberOfParticles(), positionIsSafe ? 2 : 1);
  // if the position is considered safe for insertion there should be no delayed particle
  EXPECT_EQ(delayedParticles.empty(), positionIsSafe);
}

// Generate tests for all configuration combinations
std::vector<ParameterTuple> generateParameters() {
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::mulScalar;

  // alias for readability
  const auto &particlePos = SimulationTest::testCheckedInsertParticlePos;
  const auto &searchBoxLength = SimulationTest::constellationCutoff;
  std::array<double, 3> vectorToCorner{searchBoxLength, searchBoxLength, searchBoxLength};

  return {
      // same pos as other particle
      std::make_tuple(particlePos, false),
      // pos outside iterator
      std::make_tuple(add(particlePos, {searchBoxLength * 2., 0., 0.}), true),
      // pos on sphere
      std::make_tuple(add(particlePos, {searchBoxLength, 0., 0.}), false),
      // pos outside sphere
      // from the center to the corner of a box 1/sqrt(3) = 57.74% percent of the line are within the cutoff sphere
      std::make_tuple(add(particlePos, mulScalar(vectorToCorner, 0.6)), true),
      // pos inside sphere
      std::make_tuple(add(particlePos, mulScalar(vectorToCorner, 0.4)), false),
  };
}

INSTANTIATE_TEST_SUITE_P(Generated,
                         SimulationTest,
                         testing::ValuesIn(generateParameters()),
                         SimulationTest::PrintToStringParamName());
