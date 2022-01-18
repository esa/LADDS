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
  constexpr double cutoff{1.5};
  constexpr bool newton3{false};
  constexpr size_t numDebris{3};

  std::vector<Particle> debris;
  debris.reserve(numDebris);

  // Add three particles in a row on the X axis separated by 1
  for (size_t i = 0; i < numDebris; ++i) {
    debris.emplace_back(std::array<double, 3>{static_cast<double>(i), 0., 0.},
                        std::array<double, 3>{0., static_cast<double>(i), 0.},
                        i,
                        Particle::ActivityState::passive);
  }

  CollisionFunctor collisionFunctor(cutoff, 10.0, cutoff, 0.01);

  for (size_t i = 0; i < debris.size(); ++i) {
    for (size_t j = i + 1; j < debris.size(); ++j) {
      collisionFunctor.AoSFunctor(debris[i], debris[j], newton3);
    }
  }

  // needed to merge the functor's internal thread buffers
  collisionFunctor.endTraversal(newton3);

  auto collisions = collisionFunctor.getCollisions();

  decltype(collisions) expected{
      {&debris[0], &debris[1], 1.0},
      {&debris[1], &debris[2], 1.0},
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
  debris.emplace_back(
      std::array<double, 3>{0., 0., 0.}, std::array<double, 3>{0., 0., 0.}, 0, Particle::ActivityState::passive);
  debris.back().setRadius(0.001);
  // passive small 2
  debris.emplace_back(
      std::array<double, 3>{0., 0., 0.1}, std::array<double, 3>{0., 0., 0.}, 0, Particle::ActivityState::passive);
  debris.back().setRadius(0.001);
  // passive large 1
  debris.emplace_back(
      std::array<double, 3>{0.1, 0., 0.}, std::array<double, 3>{0., 0., 0.}, 1, Particle::ActivityState::passive);
  // passive large 2
  debris.emplace_back(
      std::array<double, 3>{0.1, 0., 0.1}, std::array<double, 3>{0., 0., 0.}, 1, Particle::ActivityState::passive);
  // evasive 1
  debris.emplace_back(
      std::array<double, 3>{0., 0.1, 0.}, std::array<double, 3>{0., 0., 0.}, 2, Particle::ActivityState::evasive);
  // evasive 2
  debris.emplace_back(
      std::array<double, 3>{0., 0.1, 0.1}, std::array<double, 3>{0., 0., 0.}, 3, Particle::ActivityState::evasive);
  // evasivePreserving 1
  debris.emplace_back(std::array<double, 3>{0.1, 0.1, 0.},
                      std::array<double, 3>{0., 0., 0.},
                      4,
                      Particle::ActivityState::evasivePreserving);
  // evasivePreserving 2
  debris.emplace_back(std::array<double, 3>{0.1, 0.1, 0.1},
                      std::array<double, 3>{0., 0., 0.},
                      5,
                      Particle::ActivityState::evasivePreserving);

  CollisionFunctor collisionFunctor(cutoff, 10.0, cutoff, 0.01);

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

TEST_P(CollisionFunctorTest, LinearInterpolationTest) {
  constexpr double cutoff{80.0};
  constexpr bool newton3{false};
  constexpr size_t numDebris{2};

  const auto &[x1, x2, v1, v2, dt, expected_dist] = GetParam();

  CollisionFunctor collisionFunctor(cutoff, dt, 0.1 * cutoff, 0.1);

  std::vector<Particle> debris;
  debris.reserve(numDebris);

  // Add two particles moving in the same direction on parallel lines
  debris.emplace_back(x1, v1, 0., Particle::ActivityState::passive);
  debris.emplace_back(x2, v2, 1., Particle::ActivityState::passive);

  for (size_t i = 0; i < debris.size(); ++i) {
    for (size_t j = i + 1; j < debris.size(); ++j) {
      collisionFunctor.AoSFunctor(debris[i], debris[j], newton3);
    }
  }

  // needed to merge the functor's internal thread buffers
  collisionFunctor.endTraversal(newton3);

  auto collisions = collisionFunctor.getCollisions();

  decltype(collisions) expected{{&debris[0], &debris[1], expected_dist}};

  EXPECT_THAT(collisionFunctor.getCollisions(), ::testing::UnorderedElementsAreArray(expected));
}

// Generate tests for all configuration combinations
INSTANTIATE_TEST_SUITE_P(Generated,
                         CollisionFunctorTest,
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
