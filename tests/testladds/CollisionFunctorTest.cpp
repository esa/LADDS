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
                        i);
  }

  CollisionFunctor collisionFunctor(cutoff, 10.0, cutoff);

  for (auto &di : debris) {
    for (auto &dj : debris) {
      if (di == dj) {
        continue;
      }
      collisionFunctor.AoSFunctor(di, dj, newton3);
    }
  }

  // needed to merge the functor's internal thread buffers
  collisionFunctor.endTraversal(newton3);

  auto collisions = collisionFunctor.getCollisions();

  decltype(collisions) expected{
      {&debris[0], {&debris[1], 1.0}},
      {&debris[1], {&debris[2], 1.0}},
  };

  EXPECT_THAT(collisionFunctor.getCollisions(), ::testing::UnorderedElementsAreArray(expected));
}

TEST_P(CollisionFunctorTest, LinearInterpolationTest) {
  constexpr double cutoff{80.0};
  constexpr bool newton3{false};
  constexpr size_t numDebris{2};

  const auto &[x1, x2, v1, v2, dt, expected_dist] = GetParam();

  CollisionFunctor collisionFunctor(cutoff, dt, 0.1 * cutoff);

  std::vector<Particle> debris;
  debris.reserve(numDebris);

  // Add two particles moving in the same direction on parallel lines
  debris.emplace_back(x1, v1, 0.);
  debris.emplace_back(x2, v2, 1.);

  for (auto &di : debris) {
    for (auto &dj : debris) {
      if (di == dj) {
        continue;
      }
      collisionFunctor.AoSFunctor(di, dj, newton3);
    }
  }

  // needed to merge the functor's internal thread buffers
  collisionFunctor.endTraversal(newton3);

  auto collisions = collisionFunctor.getCollisions();

  decltype(collisions) expected{{&debris[0], {&debris[1], expected_dist}}};

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
