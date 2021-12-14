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

#include "ladds/Simulation.h"
#include "ladds/io/Logger.h"
#include "ladds/particle/Particle.h"

/**
 * Tests whether particles are correctly inserted into the simulation, when a particle
 * from the simulation and from the insertion overlap (or are very very close to each other)
 */
TEST_F(SimulationTest, testInsertionOverlap) {
  Logger logger{"SimulationTestLogger"};
  logger.get()->set_level(Logger::Level::trace);
  Simulation simulation(logger);

  auto autopas = simulation.initAutoPas(config);
  auto [csvWriter, accumulator, integrator] = simulation.initIntegrator(*autopas, config);

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

  Particle p1 = Particle({6871, 0, 0}, {0, 4.8958309146899, 5.83462408131549}, 1);
  Particle p2 = Particle({6870.99577848984, 4.89582991243564, 5.83462288687537},
                         {-0.00844301944979238, 4.89582790671436, 5.83462049654983},
                         2);
  Particle p3 = Particle({6871, 0, 0}, {0, 4.8958309146899, 5.83462408131549}, 1);
  Particle p4 = Particle({6870.99577848984, 4.89582991243564, 5.83462288687537},
                         {-0.00844301944979238, 4.89582790671436, 5.83462049654983},
                         4);
  std::vector<Particle> newSatellites;
  newSatellites.push_back(p3);
  newSatellites.push_back(p4);

  autopas->addParticle(p1);
  autopas->addParticle(p2);
  auto escapedParticles = autopas->updateContainer();
  EXPECT_EQ(escapedParticles.size(), 0);
  ASSERT_EQ(autopas->getNumberOfParticles(), 2) << "Container initialized wrong!";

  // 0
  std::vector<Particle> delayedInsertion =
      simulation.checkedInsert(*autopas, newSatellites, config["io"]["constellationCutoff"].as<double>());
  EXPECT_EQ(delayedInsertion.size(), 2);
  EXPECT_EQ(autopas->getNumberOfParticles(), 2);
  // 1
  integrator->integrate();
  escapedParticles = autopas->updateContainer();
  EXPECT_EQ(escapedParticles.size(), 0);
  delayedInsertion =
      simulation.checkedInsert(*autopas, delayedInsertion, config["io"]["constellationCutoff"].as<double>());
  EXPECT_EQ(delayedInsertion.size(), 1);
  EXPECT_EQ(autopas->getNumberOfParticles(), 3);
  // 2
  integrator->integrate();
  escapedParticles = autopas->updateContainer();
  EXPECT_EQ(escapedParticles.size(), 0);
  delayedInsertion =
      simulation.checkedInsert(*autopas, delayedInsertion, config["io"]["constellationCutoff"].as<double>());
  EXPECT_EQ(delayedInsertion.size(), 1);
  EXPECT_EQ(autopas->getNumberOfParticles(), 3);
  // 3
  integrator->integrate();
  escapedParticles = autopas->updateContainer();
  EXPECT_EQ(escapedParticles.size(), 0);
  delayedInsertion =
      simulation.checkedInsert(*autopas, delayedInsertion, config["io"]["constellationCutoff"].as<double>());
  EXPECT_EQ(delayedInsertion.size(), 0);
  EXPECT_EQ(autopas->getNumberOfParticles(), 4);
  // 4
  integrator->integrate();
  escapedParticles = autopas->updateContainer();
  EXPECT_EQ(escapedParticles.size(), 0);
  delayedInsertion =
      simulation.checkedInsert(*autopas, delayedInsertion, config["io"]["constellationCutoff"].as<double>());
  EXPECT_EQ(delayedInsertion.size(), 0);
  EXPECT_EQ(autopas->getNumberOfParticles(), 4);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 1: particle in box iterator but outside the range; is expected to be added
 */
TEST_F(SimulationTest, testCriticalRangeInsertion1) {
  Logger logger{"SimulationTestLogger"};
  logger.get()->set_level(Logger::Level::trace);
  Simulation simulation(logger);

  auto autopas = simulation.initAutoPas(config);

  std::array<double, 3> position0 = {6871, 0, 0};
  Particle p0 = Particle(position0, {0, 4.8958309146899, 5.83462408131549}, 0);

  // vector {+-length,+-length,+-length} has norm 0.04
  double length = sqrt(0.04 * 0.04 / 3.0);
  // midpoint + borderVector is on the surface of the collision sphere
  std::array<double, 3> borderVector = {length, length, length};
  std::array<double, 3> epsilon = {0.00001, 0.00001, 0.00001};

  std::array<double, 3> midPoint = position0;
  // outside : midpoint + ((l+eps,l+eps,l+eps))
  std::array<double, 3> position1 = autopas::utils::ArrayMath::add(
      midPoint, autopas::utils::ArrayMath::mul({1, 1, 1}, autopas::utils::ArrayMath::add(borderVector, epsilon)));

  Particle p1 = Particle(position1, {0, 0, 0}, 1);

  autopas->addParticle(p0);
  const auto escapedParticles = autopas->updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p1);

  std::vector<Particle> delayedInsertion =
      simulation.checkedInsert(*autopas, newSatellites, config["io"]["constellationCutoff"].as<double>());
  EXPECT_EQ(autopas->getNumberOfParticles(), 2);
  EXPECT_EQ(delayedInsertion.size(), 0);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 2: particle in box iterator and within range; is expected to be delayed
 */
TEST_F(SimulationTest, testCriticalRangeInsertion2) {
  Logger logger{"SimulationTestLogger"};
  logger.get()->set_level(Logger::Level::trace);
  Simulation simulation(logger);

  auto autopas = simulation.initAutoPas(config);

  std::array<double, 3> position0 = {6871, 0, 0};
  Particle p0 = Particle(position0, {0, 4.8958309146899, 5.83462408131549}, 0);

  // vector {+-length,+-length,+-length} has norm 0.04
  double length = sqrt(0.04 * 0.04 / 3.0);
  std::cout << length << std::endl;
  // midpoint + borderVector is on the surface of the collision sphere
  std::array<double, 3> borderVector = {length, length, length};
  std::array<double, 3> epsilon = {0.00001, 0.00001, 0.00001};

  std::array<double, 3> midPoint = position0;

  // inside : midPoint + ((l-eps,l-eps,l-eps))
  std::array<double, 3> position2 = autopas::utils::ArrayMath::add(
      midPoint, autopas::utils::ArrayMath::mul({1, 1, 1}, autopas::utils::ArrayMath::sub(borderVector, epsilon)));

  Particle p2 = Particle(position2, {0, 0, 0}, 2);

  autopas->addParticle(p0);
  const auto escapedParticles = autopas->updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p2);

  std::vector<Particle> delayedInsertion =
      simulation.checkedInsert(*autopas, newSatellites, config["io"]["constellationCutoff"].as<double>());
  EXPECT_EQ(autopas->getNumberOfParticles(), 1);
  EXPECT_EQ(delayedInsertion.size(), 1);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 3: particle in box iterator but outside of the range; is expected to be added
 */
TEST_F(SimulationTest, testCriticalRangeInsertion3) {
  Logger logger{"SimulationTestLogger"};
  logger.get()->set_level(Logger::Level::trace);
  Simulation simulation(logger);

  auto autopas = simulation.initAutoPas(config);

  std::array<double, 3> position0 = {6871, 0, 0};
  Particle p0 = Particle(position0, {0, 4.8958309146899, 5.83462408131549}, 0);

  // vector {+-length,+-length,+-length} has norm 0.04
  double length = sqrt(0.04 * 0.04 / 3.0);
  std::cout << length << std::endl;
  // midpoint + borderVector is on the surface of the collision sphere
  std::array<double, 3> borderVector = {length, length, length};
  std::array<double, 3> epsilon = {0.00001, 0.00001, 0.00001};

  std::array<double, 3> midPoint = position0;

  // outside : midPoint + ((l+eps,l+eps,l+eps)*(-1,-1,1))
  std::array<double, 3> position3 = autopas::utils::ArrayMath::add(
      midPoint, autopas::utils::ArrayMath::mul({-1, -1, 1}, autopas::utils::ArrayMath::add(borderVector, epsilon)));

  Particle p3 = Particle(position3, {0, 0, 0}, 3);

  autopas->addParticle(p0);
  const auto escapedParticles = autopas->updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p3);

  std::vector<Particle> delayedInsertion =
      simulation.checkedInsert(*autopas, newSatellites, config["io"]["constellationCutoff"].as<double>());
  EXPECT_EQ(autopas->getNumberOfParticles(), 2);
  EXPECT_EQ(delayedInsertion.size(), 0);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 4: particle in box iterator and within range; is expected to be delayed
 */
TEST_F(SimulationTest, testCriticalRangeInsertion4) {
  Logger logger{"SimulationTestLogger"};
  logger.get()->set_level(Logger::Level::trace);
  Simulation simulation(logger);

  auto autopas = simulation.initAutoPas(config);

  std::array<double, 3> position0 = {6871, 0, 0};
  Particle p0 = Particle(position0, {0, 4.8958309146899, 5.83462408131549}, 0);

  // vector {+-length,+-length,+-length} has norm 0.04
  double length = sqrt(0.04 * 0.04 / 3.0);
  std::cout << length << std::endl;
  // midpoint + borderVector is on the surface of the collision sphere
  std::array<double, 3> borderVector = {length, length, length};
  std::array<double, 3> epsilon = {0.00001, 0.00001, 0.00001};

  std::array<double, 3> midPoint = position0;

  // inside : midPoint + ((l-eps,l-eps,l-eps)*(-1,-1,1))
  std::array<double, 3> position4 = autopas::utils::ArrayMath::add(
      midPoint, autopas::utils::ArrayMath::mul({-1, -1, 1}, autopas::utils::ArrayMath::sub(borderVector, epsilon)));

  Particle p4 = Particle(position4, {0, 0, 0}, 4);

  autopas->addParticle(p0);
  const auto escapedParticles = autopas->updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p4);

  std::vector<Particle> delayedInsertion =
      simulation.checkedInsert(*autopas, newSatellites, config["io"]["constellationCutoff"].as<double>());
  EXPECT_EQ(autopas->getNumberOfParticles(), 1);
  EXPECT_EQ(delayedInsertion.size(), 1);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 5: particle in box iterator and within the range; is expected to be delayed
 */
TEST_F(SimulationTest, testCriticalRangeInsertion5) {
  Logger logger{"SimulationTestLogger"};
  logger.get()->set_level(Logger::Level::trace);
  Simulation simulation(logger);

  auto autopas = simulation.initAutoPas(config);

  std::array<double, 3> position0 = {6871, 0, 0};
  Particle p0 = Particle(position0, {0, 4.8958309146899, 5.83462408131549}, 0);

  // vector {+-length,+-length,+-length} has norm 0.04
  double length = sqrt(0.04 * 0.04 / 3.0);
  std::cout << length << std::endl;
  // midpoint + borderVector is on the surface of the collision sphere
  std::array<double, 3> borderVector = {length, length, length};
  std::array<double, 3> epsilon = {0.00001, 0.00001, 0.00001};

  std::array<double, 3> midPoint = position0;

  // inside : midPoint + ((l-eps,l-eps,l-eps)*(1,-1,1))
  std::array<double, 3> position5 = autopas::utils::ArrayMath::add(
      midPoint, autopas::utils::ArrayMath::mul({1, -1, 1}, autopas::utils::ArrayMath::sub(borderVector, epsilon)));

  Particle p5 = Particle(position5, {0, 0, 0}, 5);

  autopas->addParticle(p0);
  const auto escapedParticles = autopas->updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p5);

  std::vector<Particle> delayedInsertion =
      simulation.checkedInsert(*autopas, newSatellites, config["io"]["constellationCutoff"].as<double>());
  EXPECT_EQ(autopas->getNumberOfParticles(), 1);
  EXPECT_EQ(delayedInsertion.size(), 1);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 6: particle in box iterator and within range; is expected to be delayed
 */
TEST_F(SimulationTest, testCriticalRangeInsertion6) {
  Logger logger{"SimulationTestLogger"};
  logger.get()->set_level(Logger::Level::trace);
  Simulation simulation(logger);

  auto autopas = simulation.initAutoPas(config);

  std::array<double, 3> position0 = {6871, 0, 0};
  Particle p0 = Particle(position0, {0, 4.8958309146899, 5.83462408131549}, 0);

  // vector {+-length,+-length,+-length} has norm 0.04
  double length = sqrt(0.04 * 0.04 / 3.0);
  std::cout << length << std::endl;
  // midpoint + borderVector is on the surface of the collision sphere
  std::array<double, 3> borderVector = {length, length, length};
  std::array<double, 3> epsilon = {0.00001, 0.00001, 0.00001};

  std::array<double, 3> midPoint = position0;

  // inside : midPoint + ((l-eps,l-eps,l-eps)*(-1,-1,-1))
  std::array<double, 3> position6 = autopas::utils::ArrayMath::add(
      midPoint, autopas::utils::ArrayMath::mul({1, -1, 1}, autopas::utils::ArrayMath::sub(borderVector, epsilon)));

  Particle p6 = Particle(position6, {0, 0, 0}, 6);

  autopas->addParticle(p0);
  const auto escapedParticles = autopas->updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p6);

  std::vector<Particle> delayedInsertion =
      simulation.checkedInsert(*autopas, newSatellites, config["io"]["constellationCutoff"].as<double>());
  EXPECT_EQ(autopas->getNumberOfParticles(), 1);
  EXPECT_EQ(delayedInsertion.size(), 1);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 7: particle at the same place as the one in the simulation; is expected to be delayed
 */
TEST_F(SimulationTest, testCriticalRangeInsertion7) {
  Logger logger{"SimulationTestLogger"};
  logger.get()->set_level(Logger::Level::trace);
  Simulation simulation(logger);

  auto autopas = simulation.initAutoPas(config);

  std::array<double, 3> position0 = {6871, 0, 0};
  Particle p0 = Particle(position0, {0, 4.8958309146899, 5.83462408131549}, 0);

  // vector {+-length,+-length,+-length} has norm 0.04
  double length = sqrt(0.04 * 0.04 / 3.0);
  std::cout << length << std::endl;
  // midpoint + borderVector is on the surface of the collision sphere
  std::array<double, 3> borderVector = {length, length, length};
  std::array<double, 3> epsilon = {0.00001, 0.00001, 0.00001};

  std::array<double, 3> midPoint = position0;

  // inside: same location
  std::array<double, 3> position7 = position0;

  Particle p7 = Particle(position7, {0, 0, 0}, 7);

  autopas->addParticle(p0);
  const auto escapedParticles = autopas->updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p7);

  std::vector<Particle> delayedInsertion =
      simulation.checkedInsert(*autopas, newSatellites, config["io"]["constellationCutoff"].as<double>());
  EXPECT_EQ(autopas->getNumberOfParticles(), 1);
  EXPECT_EQ(delayedInsertion.size(), 1);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 8: particle in box iterator and within range; is expected to be delayed
 */
TEST_F(SimulationTest, testCriticalRangeInsertion8) {
  Logger logger{"SimulationTestLogger"};
  logger.get()->set_level(Logger::Level::trace);
  Simulation simulation(logger);

  auto autopas = simulation.initAutoPas(config);

  std::array<double, 3> position0 = {6871, 0, 0};
  Particle p0 = Particle(position0, {0, 4.8958309146899, 5.83462408131549}, 0);

  // vector {+-length,+-length,+-length} has norm 0.04
  double length = sqrt(0.04 * 0.04 / 3.0);
  std::cout << length << std::endl;
  // midpoint + borderVector is on the surface of the collision sphere
  std::array<double, 3> borderVector = {length, length, length};
  std::array<double, 3> epsilon = {0.00001, 0.00001, 0.00001};

  std::array<double, 3> midPoint = position0;

  // inside : midPoint + (0,0.003,0))
  std::array<double, 3> position8 = autopas::utils::ArrayMath::add(midPoint, {0, 0.003, 0});

  Particle p8 = Particle(position8, {0, 0, 0}, 8);

  autopas->addParticle(p0);
  const auto escapedParticles = autopas->updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p8);

  std::vector<Particle> delayedInsertion =
      simulation.checkedInsert(*autopas, newSatellites, config["io"]["constellationCutoff"].as<double>());
  EXPECT_EQ(autopas->getNumberOfParticles(), 1);
  EXPECT_EQ(delayedInsertion.size(), 1);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 9: particle outside box iterator and outside of the range; is expected to be added
 */
TEST_F(SimulationTest, testCriticalRangeInsertion9) {
  Logger logger{"SimulationTestLogger"};
  logger.get()->set_level(Logger::Level::trace);
  Simulation simulation(logger);

  auto autopas = simulation.initAutoPas(config);

  std::array<double, 3> position0 = {6871, 0, 0};
  Particle p0 = Particle(position0, {0, 4.8958309146899, 5.83462408131549}, 0);

  // vector {+-length,+-length,+-length} has norm 0.04
  double length = sqrt(0.04 * 0.04 / 3.0);
  std::cout << length << std::endl;
  // midpoint + borderVector is on the surface of the collision sphere
  std::array<double, 3> borderVector = {length, length, length};
  std::array<double, 3> epsilon = {0.00001, 0.00001, 0.00001};

  std::array<double, 3> midPoint = position0;

  // outside the box : midPoint + (0,0.05,0)
  std::array<double, 3> position9 = autopas::utils::ArrayMath::add(midPoint, {0, 0.05, 0});

  Particle p9 = Particle(position9, {0, 0, 0}, 9);

  autopas->addParticle(p0);
  const auto escapedParticles = autopas->updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p9);

  std::vector<Particle> delayedInsertion =
      simulation.checkedInsert(*autopas, newSatellites, config["io"]["constellationCutoff"].as<double>());
  EXPECT_EQ(autopas->getNumberOfParticles(), 2);
  EXPECT_EQ(delayedInsertion.size(), 0);
}
