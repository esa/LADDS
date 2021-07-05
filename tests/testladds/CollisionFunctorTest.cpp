/**
 * @file CollisionFunctorTest.cpp.cc
 * @author F. Gratl
 * @date 05/07/2021
 */

#include "CollisionFunctorTest.h"

#include <gmock/gmock-matchers.h>

#include "ladds/CollisionFunctor.h"
#include "ladds/Debris.h"

/**
 * Set up three particles in a row.
 * The functor is expected to identify a collision between each particle and its direct neighbors.
 */
TEST(CollisionFunctorTest, ThreeParticles) {
  constexpr double cutoff{1.5};
  constexpr size_t numDebris{3};

  std::vector<Debris> debris;
  debris.reserve(numDebris);

  // Add three particles in a row on the X axis separated by 1
  for (size_t i = 0; i < numDebris; ++i) {
    debris.emplace_back(std::array<double, 3>{static_cast<double>(i), 0., 0.}, std::array<double, 3>{0., 0., 0.}, i);
  }

  CollisionFunctor collisionFunctor(cutoff);

  for (auto &di : debris) {
    for (auto &dj : debris) {
      if (di == dj) {
        continue;
      }
      collisionFunctor.AoSFunctor(di, dj, false);
    }
  }

  auto collisions = collisionFunctor.getCollisions();

  decltype(collisions) expected{
      {&debris[0], &debris[1]},
      {&debris[1], &debris[0]},
      {&debris[1], &debris[2]},
      {&debris[2], &debris[1]},
  };

  EXPECT_THAT(collisionFunctor.getCollisions(), ::testing::UnorderedElementsAreArray(expected));
}