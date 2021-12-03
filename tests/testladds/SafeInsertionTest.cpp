//
// Created by albert on 27.11.21.
//

#include "SafeInsertionTest.h"

#include <autopas/AutoPasDecl.h>
#include <satellitePropagator/io/FileOutput.h>
#include <satellitePropagator/physics/AccelerationAccumulator.h>
#include <satellitePropagator/physics/Integrator.h>

#include "../../src/ladds/SafeInsertion.h"

/**
 * Tests whether particles are correctly inserted into the simulation, when a particle
 * from the simulation and from the insertion overlap (or are very very close to each other)
 */
TEST(SafeInsertionTest, testInsertionOverlap) {
  using AutoPas_t = autopas::AutoPas<Particle>;

  autopas::AutoPas<Particle> autopas;
  std::array<bool, 8> selectedPropagatorComponents = {true, false, false, false, false, false, false, false};

  auto fo =
      std::make_shared<FileOutput<AutoPas_t>>(autopas, "output.csv", OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t>>(
      selectedPropagatorComponents, autopas, 0.0, *fo);
  auto deltaT = 1.0;
  auto integrator = std::make_shared<Integrator<AutoPas_t>>(autopas, *accumulator, deltaT);

  const auto maxAltitude = 85000;
  const auto cutoff = 0.02;
  const auto verletSkin = 0.20;
  const auto verletRebuildFrequency = 20;
  const auto desiredCellsPerDimension = 32;

  autopas.setBoxMin({-maxAltitude, -maxAltitude, -maxAltitude});
  autopas.setBoxMax({maxAltitude, maxAltitude, maxAltitude});
  autopas.setCutoff(cutoff);
  autopas.setVerletSkin(verletSkin);
  autopas.setVerletRebuildFrequency(verletRebuildFrequency);
  autopas.setCellSizeFactor((maxAltitude * 2.) / ((cutoff + verletSkin) * (desiredCellsPerDimension - 2)));
  autopas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
  autopas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autopas.setAllowedContainers({autopas::ContainerOption::varVerletListsAsBuild});
  autopas.setAllowedTraversals({autopas::TraversalOption::vvl_as_built});
  autopas::Logger::get()->set_level(spdlog::level::from_str("off"));
  autopas.init();

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

  autopas.addParticle(p1);
  autopas.addParticle(p2);
  const auto escapedParticles = autopas.updateContainer();

  // 0
  std::vector<Particle> delayedInsertion = SafeInsertion::insert(autopas, newSatellites, cutoff);
  ASSERT_EQ(delayedInsertion.size(), 2);
  ASSERT_EQ(autopas.getNumberOfParticles(), 2);
  // 1
  integrator->integrate();
  delayedInsertion = SafeInsertion::insert(autopas, delayedInsertion, cutoff);
  ASSERT_EQ(delayedInsertion.size(), 1);
  ASSERT_EQ(autopas.getNumberOfParticles(), 3);
  // 2
  integrator->integrate();
  delayedInsertion = SafeInsertion::insert(autopas, delayedInsertion, cutoff);
  ASSERT_EQ(delayedInsertion.size(), 1);
  ASSERT_EQ(autopas.getNumberOfParticles(), 3);
  // 3
  integrator->integrate();
  delayedInsertion = SafeInsertion::insert(autopas, delayedInsertion, cutoff);
  ASSERT_EQ(delayedInsertion.size(), 0);
  ASSERT_EQ(autopas.getNumberOfParticles(), 4);
  // 4
  integrator->integrate();
  delayedInsertion = SafeInsertion::insert(autopas, delayedInsertion, cutoff);
  ASSERT_EQ(delayedInsertion.size(), 0);
  ASSERT_EQ(autopas.getNumberOfParticles(), 4);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 1: particle in box iterator but outside the range; is expected to be added
 */
TEST(SafeInsertionTest, testCriticalRangeInsertion1) {
  using AutoPas_t = autopas::AutoPas<Particle>;

  autopas::AutoPas<Particle> autopas;
  std::array<bool, 8> selectedPropagatorComponents = {true, false, false, false, false, false, false, false};

  auto fo =
      std::make_shared<FileOutput<AutoPas_t>>(autopas, "output.csv", OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t>>(
      selectedPropagatorComponents, autopas, 0.0, *fo);
  auto deltaT = 1.0;
  auto integrator = std::make_shared<Integrator<AutoPas_t>>(autopas, *accumulator, deltaT);

  const auto maxAltitude = 85000;
  const auto cutoff = 0.02;
  const auto verletSkin = 0.20;
  const auto verletRebuildFrequency = 20;
  const auto desiredCellsPerDimension = 32;

  autopas.setBoxMin({-maxAltitude, -maxAltitude, -maxAltitude});
  autopas.setBoxMax({maxAltitude, maxAltitude, maxAltitude});
  autopas.setCutoff(cutoff);
  autopas.setVerletSkin(verletSkin);
  autopas.setVerletRebuildFrequency(verletRebuildFrequency);
  autopas.setCellSizeFactor((maxAltitude * 2.) / ((cutoff + verletSkin) * (desiredCellsPerDimension - 2)));
  autopas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
  autopas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autopas.setAllowedContainers({autopas::ContainerOption::varVerletListsAsBuild});
  autopas.setAllowedTraversals({autopas::TraversalOption::vvl_as_built});
  autopas::Logger::get()->set_level(spdlog::level::from_str("off"));
  autopas.init();

  std::array<double, 3> position0 = {6871, 0, 0};
  Particle p0 = Particle(position0, {0, 4.8958309146899, 5.83462408131549}, 0);

  // vector {+-length,+-length,+-length} has norm 0.04
  double length = sqrt(0.04 * 0.04 / 3.0);
  std::cout << length << std::endl;
  // midpoint + borderVector is on the surface of the collision sphere
  std::array<double, 3> borderVector = {length, length, length};
  std::array<double, 3> epsilon = {0.00001, 0.00001, 0.00001};

  std::array<double, 3> midPoint = position0;
  // outside : midpoint + ((l+eps,l+eps,l+eps))
  std::array<double, 3> position1 = autopas::utils::ArrayMath::add(
      midPoint, autopas::utils::ArrayMath::mul({1, 1, 1}, autopas::utils::ArrayMath::add(borderVector, epsilon)));

  Particle p1 = Particle(position1, {0, 0, 0}, 1);

  autopas.addParticle(p0);
  const auto escapedParticles = autopas.updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p1);

  std::vector<Particle> delayedInsertion = SafeInsertion::insert(autopas, newSatellites, cutoff);
  ASSERT_EQ(autopas.getNumberOfParticles(), 2);
  ASSERT_EQ(delayedInsertion.size(), 0);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 2: particle in box iterator and within range; is expected to be delayed
 */
TEST(SafeInsertionTest, testCriticalRangeInsertion2) {
  using AutoPas_t = autopas::AutoPas<Particle>;

  autopas::AutoPas<Particle> autopas;
  std::array<bool, 8> selectedPropagatorComponents = {true, false, false, false, false, false, false, false};

  auto fo =
      std::make_shared<FileOutput<AutoPas_t>>(autopas, "output.csv", OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t>>(
      selectedPropagatorComponents, autopas, 0.0, *fo);
  auto deltaT = 1.0;
  auto integrator = std::make_shared<Integrator<AutoPas_t>>(autopas, *accumulator, deltaT);

  const auto maxAltitude = 85000;
  const auto cutoff = 0.02;
  const auto verletSkin = 0.20;
  const auto verletRebuildFrequency = 20;
  const auto desiredCellsPerDimension = 32;

  autopas.setBoxMin({-maxAltitude, -maxAltitude, -maxAltitude});
  autopas.setBoxMax({maxAltitude, maxAltitude, maxAltitude});
  autopas.setCutoff(cutoff);
  autopas.setVerletSkin(verletSkin);
  autopas.setVerletRebuildFrequency(verletRebuildFrequency);
  autopas.setCellSizeFactor((maxAltitude * 2.) / ((cutoff + verletSkin) * (desiredCellsPerDimension - 2)));
  autopas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
  autopas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autopas.setAllowedContainers({autopas::ContainerOption::varVerletListsAsBuild});
  autopas.setAllowedTraversals({autopas::TraversalOption::vvl_as_built});
  autopas::Logger::get()->set_level(spdlog::level::from_str("off"));
  autopas.init();

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

  autopas.addParticle(p0);
  const auto escapedParticles = autopas.updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p2);

  std::vector<Particle> delayedInsertion = SafeInsertion::insert(autopas, newSatellites, cutoff);
  ASSERT_EQ(autopas.getNumberOfParticles(), 1);
  ASSERT_EQ(delayedInsertion.size(), 1);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 3: particle in box iterator but outside of the range; is expected to be added
 */
TEST(SafeInsertionTest, testCriticalRangeInsertion3) {
  using AutoPas_t = autopas::AutoPas<Particle>;

  autopas::AutoPas<Particle> autopas;
  std::array<bool, 8> selectedPropagatorComponents = {true, false, false, false, false, false, false, false};

  auto fo =
      std::make_shared<FileOutput<AutoPas_t>>(autopas, "output.csv", OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t>>(
      selectedPropagatorComponents, autopas, 0.0, *fo);
  auto deltaT = 1.0;
  auto integrator = std::make_shared<Integrator<AutoPas_t>>(autopas, *accumulator, deltaT);

  const auto maxAltitude = 85000;
  const auto cutoff = 0.02;
  const auto verletSkin = 0.20;
  const auto verletRebuildFrequency = 20;
  const auto desiredCellsPerDimension = 32;

  autopas.setBoxMin({-maxAltitude, -maxAltitude, -maxAltitude});
  autopas.setBoxMax({maxAltitude, maxAltitude, maxAltitude});
  autopas.setCutoff(cutoff);
  autopas.setVerletSkin(verletSkin);
  autopas.setVerletRebuildFrequency(verletRebuildFrequency);
  autopas.setCellSizeFactor((maxAltitude * 2.) / ((cutoff + verletSkin) * (desiredCellsPerDimension - 2)));
  autopas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
  autopas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autopas.setAllowedContainers({autopas::ContainerOption::varVerletListsAsBuild});
  autopas.setAllowedTraversals({autopas::TraversalOption::vvl_as_built});
  autopas::Logger::get()->set_level(spdlog::level::from_str("off"));
  autopas.init();

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

  autopas.addParticle(p0);
  const auto escapedParticles = autopas.updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p3);

  std::vector<Particle> delayedInsertion = SafeInsertion::insert(autopas, newSatellites, cutoff);
  ASSERT_EQ(autopas.getNumberOfParticles(), 2);
  ASSERT_EQ(delayedInsertion.size(), 0);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 4: particle in box iterator and within range; is expected to be delayed
 */
TEST(SafeInsertionTest, testCriticalRangeInsertion4) {
  using AutoPas_t = autopas::AutoPas<Particle>;

  autopas::AutoPas<Particle> autopas;
  std::array<bool, 8> selectedPropagatorComponents = {true, false, false, false, false, false, false, false};

  auto fo =
      std::make_shared<FileOutput<AutoPas_t>>(autopas, "output.csv", OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t>>(
      selectedPropagatorComponents, autopas, 0.0, *fo);
  auto deltaT = 1.0;
  auto integrator = std::make_shared<Integrator<AutoPas_t>>(autopas, *accumulator, deltaT);

  const auto maxAltitude = 85000;
  const auto cutoff = 0.02;
  const auto verletSkin = 0.20;
  const auto verletRebuildFrequency = 20;
  const auto desiredCellsPerDimension = 32;

  autopas.setBoxMin({-maxAltitude, -maxAltitude, -maxAltitude});
  autopas.setBoxMax({maxAltitude, maxAltitude, maxAltitude});
  autopas.setCutoff(cutoff);
  autopas.setVerletSkin(verletSkin);
  autopas.setVerletRebuildFrequency(verletRebuildFrequency);
  autopas.setCellSizeFactor((maxAltitude * 2.) / ((cutoff + verletSkin) * (desiredCellsPerDimension - 2)));
  autopas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
  autopas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autopas.setAllowedContainers({autopas::ContainerOption::varVerletListsAsBuild});
  autopas.setAllowedTraversals({autopas::TraversalOption::vvl_as_built});
  autopas::Logger::get()->set_level(spdlog::level::from_str("off"));
  autopas.init();

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

  autopas.addParticle(p0);
  const auto escapedParticles = autopas.updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p4);

  std::vector<Particle> delayedInsertion = SafeInsertion::insert(autopas, newSatellites, cutoff);
  ASSERT_EQ(autopas.getNumberOfParticles(), 1);
  ASSERT_EQ(delayedInsertion.size(), 1);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 5: particle in box iterator and within the range; is expected to be delayed
 */
TEST(SafeInsertionTest, testCriticalRangeInsertion5) {
  using AutoPas_t = autopas::AutoPas<Particle>;

  autopas::AutoPas<Particle> autopas;
  std::array<bool, 8> selectedPropagatorComponents = {true, false, false, false, false, false, false, false};

  auto fo =
      std::make_shared<FileOutput<AutoPas_t>>(autopas, "output.csv", OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t>>(
      selectedPropagatorComponents, autopas, 0.0, *fo);
  auto deltaT = 1.0;
  auto integrator = std::make_shared<Integrator<AutoPas_t>>(autopas, *accumulator, deltaT);

  const auto maxAltitude = 85000;
  const auto cutoff = 0.02;
  const auto verletSkin = 0.20;
  const auto verletRebuildFrequency = 20;
  const auto desiredCellsPerDimension = 32;

  autopas.setBoxMin({-maxAltitude, -maxAltitude, -maxAltitude});
  autopas.setBoxMax({maxAltitude, maxAltitude, maxAltitude});
  autopas.setCutoff(cutoff);
  autopas.setVerletSkin(verletSkin);
  autopas.setVerletRebuildFrequency(verletRebuildFrequency);
  autopas.setCellSizeFactor((maxAltitude * 2.) / ((cutoff + verletSkin) * (desiredCellsPerDimension - 2)));
  autopas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
  autopas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autopas.setAllowedContainers({autopas::ContainerOption::varVerletListsAsBuild});
  autopas.setAllowedTraversals({autopas::TraversalOption::vvl_as_built});
  autopas::Logger::get()->set_level(spdlog::level::from_str("off"));
  autopas.init();

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

  autopas.addParticle(p0);
  const auto escapedParticles = autopas.updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p5);

  std::vector<Particle> delayedInsertion = SafeInsertion::insert(autopas, newSatellites, cutoff);
  ASSERT_EQ(autopas.getNumberOfParticles(), 1);
  ASSERT_EQ(delayedInsertion.size(), 1);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 6: particle in box iterator and within range; is expected to be delayed
 */
TEST(SafeInsertionTest, testCriticalRangeInsertion6) {
  using AutoPas_t = autopas::AutoPas<Particle>;

  autopas::AutoPas<Particle> autopas;
  std::array<bool, 8> selectedPropagatorComponents = {true, false, false, false, false, false, false, false};

  auto fo =
      std::make_shared<FileOutput<AutoPas_t>>(autopas, "output.csv", OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t>>(
      selectedPropagatorComponents, autopas, 0.0, *fo);
  auto deltaT = 1.0;
  auto integrator = std::make_shared<Integrator<AutoPas_t>>(autopas, *accumulator, deltaT);

  const auto maxAltitude = 85000;
  const auto cutoff = 0.02;
  const auto verletSkin = 0.20;
  const auto verletRebuildFrequency = 20;
  const auto desiredCellsPerDimension = 32;

  autopas.setBoxMin({-maxAltitude, -maxAltitude, -maxAltitude});
  autopas.setBoxMax({maxAltitude, maxAltitude, maxAltitude});
  autopas.setCutoff(cutoff);
  autopas.setVerletSkin(verletSkin);
  autopas.setVerletRebuildFrequency(verletRebuildFrequency);
  autopas.setCellSizeFactor((maxAltitude * 2.) / ((cutoff + verletSkin) * (desiredCellsPerDimension - 2)));
  autopas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
  autopas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autopas.setAllowedContainers({autopas::ContainerOption::varVerletListsAsBuild});
  autopas.setAllowedTraversals({autopas::TraversalOption::vvl_as_built});
  autopas::Logger::get()->set_level(spdlog::level::from_str("off"));
  autopas.init();

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

  autopas.addParticle(p0);
  const auto escapedParticles = autopas.updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p6);

  std::vector<Particle> delayedInsertion = SafeInsertion::insert(autopas, newSatellites, cutoff);
  ASSERT_EQ(autopas.getNumberOfParticles(), 1);
  ASSERT_EQ(delayedInsertion.size(), 1);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 7: particle at the same place as the one in the simulation; is expected to be delayed
 */
TEST(SafeInsertionTest, testCriticalRangeInsertion7) {
  using AutoPas_t = autopas::AutoPas<Particle>;

  autopas::AutoPas<Particle> autopas;
  std::array<bool, 8> selectedPropagatorComponents = {true, false, false, false, false, false, false, false};

  auto fo =
      std::make_shared<FileOutput<AutoPas_t>>(autopas, "output.csv", OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t>>(
      selectedPropagatorComponents, autopas, 0.0, *fo);
  auto deltaT = 1.0;
  auto integrator = std::make_shared<Integrator<AutoPas_t>>(autopas, *accumulator, deltaT);

  const auto maxAltitude = 85000;
  const auto cutoff = 0.02;
  const auto verletSkin = 0.20;
  const auto verletRebuildFrequency = 20;
  const auto desiredCellsPerDimension = 32;

  autopas.setBoxMin({-maxAltitude, -maxAltitude, -maxAltitude});
  autopas.setBoxMax({maxAltitude, maxAltitude, maxAltitude});
  autopas.setCutoff(cutoff);
  autopas.setVerletSkin(verletSkin);
  autopas.setVerletRebuildFrequency(verletRebuildFrequency);
  autopas.setCellSizeFactor((maxAltitude * 2.) / ((cutoff + verletSkin) * (desiredCellsPerDimension - 2)));
  autopas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
  autopas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autopas.setAllowedContainers({autopas::ContainerOption::varVerletListsAsBuild});
  autopas.setAllowedTraversals({autopas::TraversalOption::vvl_as_built});
  autopas::Logger::get()->set_level(spdlog::level::from_str("off"));
  autopas.init();

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

  autopas.addParticle(p0);
  const auto escapedParticles = autopas.updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p7);

  std::vector<Particle> delayedInsertion = SafeInsertion::insert(autopas, newSatellites, cutoff);
  ASSERT_EQ(autopas.getNumberOfParticles(), 1);
  ASSERT_EQ(delayedInsertion.size(), 1);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 8: particle in box iterator and within range; is expected to be delayed
 */
TEST(SafeInsertionTest, testCriticalRangeInsertion8) {
  using AutoPas_t = autopas::AutoPas<Particle>;

  autopas::AutoPas<Particle> autopas;
  std::array<bool, 8> selectedPropagatorComponents = {true, false, false, false, false, false, false, false};

  auto fo =
      std::make_shared<FileOutput<AutoPas_t>>(autopas, "output.csv", OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t>>(
      selectedPropagatorComponents, autopas, 0.0, *fo);
  auto deltaT = 1.0;
  auto integrator = std::make_shared<Integrator<AutoPas_t>>(autopas, *accumulator, deltaT);

  const auto maxAltitude = 85000;
  const auto cutoff = 0.02;
  const auto verletSkin = 0.20;
  const auto verletRebuildFrequency = 20;
  const auto desiredCellsPerDimension = 32;

  autopas.setBoxMin({-maxAltitude, -maxAltitude, -maxAltitude});
  autopas.setBoxMax({maxAltitude, maxAltitude, maxAltitude});
  autopas.setCutoff(cutoff);
  autopas.setVerletSkin(verletSkin);
  autopas.setVerletRebuildFrequency(verletRebuildFrequency);
  autopas.setCellSizeFactor((maxAltitude * 2.) / ((cutoff + verletSkin) * (desiredCellsPerDimension - 2)));
  autopas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
  autopas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autopas.setAllowedContainers({autopas::ContainerOption::varVerletListsAsBuild});
  autopas.setAllowedTraversals({autopas::TraversalOption::vvl_as_built});
  autopas::Logger::get()->set_level(spdlog::level::from_str("off"));
  autopas.init();

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

  autopas.addParticle(p0);
  const auto escapedParticles = autopas.updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p8);

  std::vector<Particle> delayedInsertion = SafeInsertion::insert(autopas, newSatellites, cutoff);
  ASSERT_EQ(autopas.getNumberOfParticles(), 1);
  ASSERT_EQ(delayedInsertion.size(), 1);
}

/**
 * tests whether particles are only inserted to the simulation, when they are outside a critical
 * range and delayed when within that critical range
 *
 * 9: particle outside box iterator and outside of the range; is expected to be added
 */
TEST(SafeInsertionTest, testCriticalRangeInsertion9) {
  using AutoPas_t = autopas::AutoPas<Particle>;

  autopas::AutoPas<Particle> autopas;
  std::array<bool, 8> selectedPropagatorComponents = {true, false, false, false, false, false, false, false};

  auto fo =
      std::make_shared<FileOutput<AutoPas_t>>(autopas, "output.csv", OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t>>(
      selectedPropagatorComponents, autopas, 0.0, *fo);
  auto deltaT = 1.0;
  auto integrator = std::make_shared<Integrator<AutoPas_t>>(autopas, *accumulator, deltaT);

  const auto maxAltitude = 85000;
  const auto cutoff = 0.02;
  const auto verletSkin = 0.20;
  const auto verletRebuildFrequency = 20;
  const auto desiredCellsPerDimension = 32;

  autopas.setBoxMin({-maxAltitude, -maxAltitude, -maxAltitude});
  autopas.setBoxMax({maxAltitude, maxAltitude, maxAltitude});
  autopas.setCutoff(cutoff);
  autopas.setVerletSkin(verletSkin);
  autopas.setVerletRebuildFrequency(verletRebuildFrequency);
  autopas.setCellSizeFactor((maxAltitude * 2.) / ((cutoff + verletSkin) * (desiredCellsPerDimension - 2)));
  autopas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
  autopas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autopas.setAllowedContainers({autopas::ContainerOption::varVerletListsAsBuild});
  autopas.setAllowedTraversals({autopas::TraversalOption::vvl_as_built});
  autopas::Logger::get()->set_level(spdlog::level::from_str("off"));
  autopas.init();

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

  autopas.addParticle(p0);
  const auto escapedParticles = autopas.updateContainer();

  std::vector<Particle> newSatellites;
  newSatellites.push_back(p9);

  std::vector<Particle> delayedInsertion = SafeInsertion::insert(autopas, newSatellites, cutoff);
  ASSERT_EQ(autopas.getNumberOfParticles(), 2);
  ASSERT_EQ(delayedInsertion.size(), 0);
}
