/**
 * @file CollisionFunctorTest.cpp.cc
 * @author F. Gratl
 * @date 05/07/2021
 */

#include "CollisionFunctorTest.h"

#include <gmock/gmock-matchers.h>

#include "ladds/CollisionFunctor.h"
#include "ladds/particle/Particle.h"

/**
 * Set up three particles in a row.
 * The functor is expected to identify a collision between each particle and its direct neighbors.
 */
TEST(CollisionFunctorTest, ThreeParticles) {
  constexpr double cutoff{1};
  constexpr bool newton3{false};
  constexpr size_t numDebris{3};

  std::vector<Particle> debris;
  debris.reserve(numDebris);

  // Add three particles with radius .6m in a row on the X axis separated by 1 m
  // => direct neighbors overlap
  for (size_t i = 0; i < numDebris; ++i) {
    debris.emplace_back(std::array<double, 3>{static_cast<double>(i) / 1000., 0., 0.},
                        std::array<double, 3>{0., static_cast<double>(i), 0.},
                        i,
                        "dummy",
                        Particle::ActivityState::passive,
                        1.,
                        .6,
                        Particle::calculateBcInv(0., 1., 1., 2.2));
  }

  CollisionFunctor collisionFunctor(cutoff, 10.0, 1., 0.01);

  for (size_t i = 0; i < debris.size(); ++i) {
    for (size_t j = i + 1; j < debris.size(); ++j) {
      collisionFunctor.AoSFunctor(debris[i], debris[j], newton3);
    }
  }

  // needed to merge the functor's internal thread buffers
  collisionFunctor.endTraversal(newton3);

  auto collisions = collisionFunctor.getCollisions();

  decltype(collisions) expected{
      {&debris[0], &debris[1], 1.0 / (1000 * 1000)},  // convert distance to km^2
      {&debris[1], &debris[2], 1.0 / (1000 * 1000)},  // convert distance to km^2
  };

  EXPECT_THAT(collisionFunctor.getCollisions(), ::testing::UnorderedElementsAreArray(expected));
}

/**
 * Set up two particles per ActivityState, all very close to each other.
 * We expect exactly the following collisions:
 *   passiveSmall      x  passiveSmall
 *   passiveSmall      x  passiveLarge
 *   passiveSmall      x  evasive
 *   passiveSmall      x  evasivePreserving
 *   passiveLarge      x  passiveLarge
 * Thus the following collisions to be disregarded:
 *   evasive           x  passiveLarge
 *   evasive           x  evasive
 *   evasive           x  evadingPreserving
 *   evasivePreserving x  passiveLarge
 *   evasivePreserving x  evadingPreserving
 */
TEST(CollisionFunctorTest, MixActivityStates) {
  constexpr double cutoff{1.};
  constexpr bool newton3{false};
  constexpr size_t numDebris{8};

  std::vector<Particle> debris;
  debris.reserve(numDebris);

  // passive small 1
  debris.emplace_back(std::array<double, 3>{0., 0., 0.},
                      std::array<double, 3>{0., 0., 0.},
                      0,
                      "passive small 1",
                      Particle::ActivityState::passive,
                      1.,
                      0.001,
                      Particle::calculateBcInv(0., 1., 1., 2.2));
  // passive small 2
  debris.emplace_back(std::array<double, 3>{0., 0., 0.1},
                      std::array<double, 3>{0., 0., 0.},
                      0,
                      "passive small 2",
                      Particle::ActivityState::passive,
                      1.,
                      0.001,
                      Particle::calculateBcInv(0., 1., 1., 2.2));
  // passive large 1
  debris.emplace_back(std::array<double, 3>{0.1, 0., 0.},
                      std::array<double, 3>{0., 0., 0.},
                      1,
                      "passive large 1",
                      Particle::ActivityState::passive,
                      1.,
                      1.,
                      Particle::calculateBcInv(0., 1., 1., 2.2));
  // passive large 2
  debris.emplace_back(std::array<double, 3>{0.1, 0., 0.1},
                      std::array<double, 3>{0., 0., 0.},
                      1,
                      "passive large 2",
                      Particle::ActivityState::passive,
                      1.,
                      1.,
                      Particle::calculateBcInv(0., 1., 1., 2.2));
  // evasive 1
  debris.emplace_back(std::array<double, 3>{0., 0.1, 0.},
                      std::array<double, 3>{0., 0., 0.},
                      2,
                      "evasive 1",
                      Particle::ActivityState::evasive,
                      1.,
                      1.,
                      Particle::calculateBcInv(0., 1., 1., 2.2));
  // evasive 2
  debris.emplace_back(std::array<double, 3>{0., 0.1, 0.1},
                      std::array<double, 3>{0., 0., 0.},
                      3,
                      "evasive 2",
                      Particle::ActivityState::evasive,
                      1.,
                      1.,
                      Particle::calculateBcInv(0., 1., 1., 2.2));
  // evasivePreserving 1
  debris.emplace_back(std::array<double, 3>{0.1, 0.1, 0.},
                      std::array<double, 3>{0., 0., 0.},
                      4,
                      "evasivePreserving 1",
                      Particle::ActivityState::evasivePreserving,
                      1.,
                      1.,
                      Particle::calculateBcInv(0., 1., 1., 2.2));
  // evasivePreserving 2
  debris.emplace_back(std::array<double, 3>{0.1, 0.1, 0.1},
                      std::array<double, 3>{0., 0., 0.},
                      5,
                      "evasivePreserving 2",
                      Particle::ActivityState::evasivePreserving,
                      1.,
                      1.,
                      Particle::calculateBcInv(0., 1., 1., 2.2));

  CollisionFunctor collisionFunctor(cutoff, 10.0, 1., 0.01);

  for (size_t i = 0; i < debris.size(); ++i) {
    for (size_t j = i + 1; j < debris.size(); ++j) {
      collisionFunctor.AoSFunctor(debris[i], debris[j], newton3);
    }
  }

  // needed to merge the functor's internal thread buffers
  collisionFunctor.endTraversal(newton3);

  auto collisions = collisionFunctor.getCollisions();

  // extract ID pairs from collisions
  std::vector<std::tuple<size_t, size_t>> collisionIdPairs;
  collisionIdPairs.reserve(collisions.size());
  std::transform(collisions.begin(), collisions.end(), std::back_inserter(collisionIdPairs), [](const auto &collision) {
    const auto &[d1, d2, dist] = collision;
    return std::make_tuple(d1->getID(), d2->getID());
  });

  // set up expectations
  decltype(collisionIdPairs) expectations{
      // passive small 1 collides with everything
      {debris[0].getID(), debris[1].getID()},
      {debris[0].getID(), debris[2].getID()},
      {debris[0].getID(), debris[3].getID()},
      {debris[0].getID(), debris[4].getID()},
      {debris[0].getID(), debris[5].getID()},
      {debris[0].getID(), debris[6].getID()},
      {debris[0].getID(), debris[7].getID()},
      // passive small 2 collides with everything (collsion with passive small 1 already covered above)
      {debris[1].getID(), debris[2].getID()},
      {debris[1].getID(), debris[3].getID()},
      {debris[1].getID(), debris[4].getID()},
      {debris[1].getID(), debris[5].getID()},
      {debris[1].getID(), debris[6].getID()},
      {debris[1].getID(), debris[7].getID()},
      // passive large 1 and 2 collide
      {debris[2].getID(), debris[3].getID()},
  };

  EXPECT_EQ(collisions.size(), expectations.size());
  EXPECT_THAT(collisionIdPairs, ::testing::UnorderedElementsAreArray(expectations));
}

/**
 * Two particles each with radius 1 are placed at a distance of 3 m.
 * Using collisionDistanceFactor 1 they should not be considered colliding, using 2 they should.
 */
TEST(CollisionFunctorTest, CollisionDistanceFactorTest) {
  constexpr double cutoff{80.0};
  constexpr double dt{1.};
  constexpr bool newton3{false};
  constexpr double particleMass{1.};
  constexpr double particleRadius{1.};

  std::vector<Particle> debris{{{0., 0., 0.},
                                {1., 0., 0.},
                                1,
                                "A",
                                Particle::ActivityState::passive,
                                particleMass,
                                particleRadius,
                                Particle::calculateBcInv(0., particleMass, particleRadius, 2.2)},
                               {{3. / 1000, 0., 0.},
                                {-1., 0., 0.},
                                2,
                                "B",
                                Particle::ActivityState::passive,
                                particleMass,
                                particleRadius,
                                Particle::calculateBcInv(0., particleMass, particleRadius, 2.2)}};

  // test for different factors
  for (const double collisionDistanceFactor : {1., 2.}) {
    CollisionFunctor collisionFunctor(cutoff, dt, collisionDistanceFactor, 0.1);
    collisionFunctor.AoSFunctor(debris[0], debris[1], newton3);
    collisionFunctor.endTraversal(newton3);
    if (collisionDistanceFactor == 1.) {
      // distance too far -> no collisions
      EXPECT_EQ(collisionFunctor.getCollisions().size(), 0)
          << "For collisionDistanceFactor: " << collisionDistanceFactor;
    } else if (collisionDistanceFactor == 2.) {
      // over approximation should now include particles -> collision
      EXPECT_EQ(collisionFunctor.getCollisions().size(), 1)
          << "For collisionDistanceFactor: " << collisionDistanceFactor;
    } else {
      GTEST_FAIL() << "Unexpected collisionDistanceFactor: " << collisionDistanceFactor;
    }
  }
}

TEST_P(CollisionFunctorTest, LinearInterpolationTest) {
  constexpr double cutoff{80.0};
  constexpr bool newton3{false};
  constexpr size_t numDebris{2};
  constexpr double particleRadius{1000.};

  const auto &[x1, x2, v1, v2, dt, squaredExpectedDist] = GetParam();

  // collisionDistanceFactor > 1 to actually cover all conjunctions
  CollisionFunctor collisionFunctor(cutoff, dt, 10., 0.1);

  std::vector<Particle> debris;
  debris.reserve(numDebris);

  // Add two particles moving in the same direction on parallel lines
  debris.emplace_back(
      x1, v1, 0., "A", Particle::ActivityState::passive, 1., particleRadius, Particle::calculateBcInv(0., 1., 1., 2.2));
  debris.emplace_back(
      x2, v2, 1., "B", Particle::ActivityState::passive, 1., particleRadius, Particle::calculateBcInv(0., 1., 1., 2.2));

  for (size_t i = 0; i < debris.size(); ++i) {
    for (size_t j = i + 1; j < debris.size(); ++j) {
      collisionFunctor.AoSFunctor(debris[i], debris[j], newton3);
    }
  }

  // needed to merge the functor's internal thread buffers
  collisionFunctor.endTraversal(newton3);

  auto collisions = collisionFunctor.getCollisions();

  decltype(collisions) expected{{&debris[0], &debris[1], squaredExpectedDist}};

  // helper function for debugging output
  auto getIDsStringFromPointers = [](const auto &collisions) {
    std::stringstream ss;
    std::vector<std::tuple<size_t, size_t>> ids;
    for (const auto &[ptrA, ptrB, dist] : collisions) {
      ss << "[" << ptrA->getID() << "|" << ptrB->getID() << "|" << dist << "]";
    }
    return ss.str();
  };

  EXPECT_THAT(collisionFunctor.getCollisions(), ::testing::UnorderedElementsAreArray(expected))
      << "Expected tuples: " << getIDsStringFromPointers(expected)
      << "\nFound    tuples: " << getIDsStringFromPointers(collisions);
}

// Generate tests for all configuration combinations
INSTANTIATE_TEST_SUITE_P(Generated,
                         CollisionFunctorTest,
                         // x1 [km], x2 [km], v1 [km/s], v2 [km/2], dt [s], squared expected_dist [km^2]
                         testing::ValuesIn({std::make_tuple(std::array<double, 3>{0., 1., 2.},
                                                            std::array<double, 3>{2., 4., 6.},
                                                            std::array<double, 3>{0., 0., 0.},
                                                            std::array<double, 3>{1., 2., 3.},
                                                            1.,
                                                            3.),
                                            std::make_tuple(std::array<double, 3>{5., 0., 0.},
                                                            std::array<double, 3>{1., 1., 1.},
                                                            std::array<double, 3>{-1., -1., -1.},
                                                            std::array<double, 3>{1., 1., 1.},
                                                            1.,
                                                            18.),
                                            std::make_tuple(std::array<double, 3>{1.5, 0., 1.},
                                                            std::array<double, 3>{3., 0., 2.},
                                                            std::array<double, 3>{-3., 0., -2.},
                                                            std::array<double, 3>{3., 0., 2.},
                                                            2.,
                                                            0.),
                                            std::make_tuple(std::array<double, 3>{0., 0., 0.},
                                                            std::array<double, 3>{0., 0.5, 0},
                                                            std::array<double, 3>{1., 0, 0},
                                                            std::array<double, 3>{0., 0, 1.},
                                                            42.,
                                                            .25)}),
                         CollisionFunctorTest::PrintToStringParamName());
