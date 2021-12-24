/**
 * @file CollisionFunctorTest.h
 * @author F. Gratl
 * @date 05/07/2021
 */

#include "CollisionFunctorIntegrationTest.h"

#include <autopas/AutoPasDecl.h>
#include <autopasTools/generators/RandomGenerator.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>
#include <ladds/CollisionFunctor.h>

extern template class autopas::AutoPas<Particle>;
extern template bool autopas::AutoPas<Particle>::iteratePairwise(CollisionFunctor *);

void CollisionFunctorIntegrationTest::SetUpTestSuite() {
  using autopasTools::generators::RandomGenerator;
  constexpr size_t numDebris = 15;
  Particle defaultParticle{{
                               0.,
                               0.,
                               0.,
                           },
                           {
                               0.,
                               0.,
                               0.,
                           },
                           0};
  _debris.reserve(numDebris);

  // fix seed for randomPosition()
  srand(42);

  for (int i = 0; i < numDebris; ++i) {
    defaultParticle.setR(autopasTools::generators::RandomGenerator::randomPosition(_boxMin, _boxMax));
    defaultParticle.setID(i);
    _debris.push_back(defaultParticle);
  }

  // if distance is smaller than cutoff save the pair
  for (auto &i : _debris) {
    for (auto &j : _debris) {
      // compare pointer since we are looking at the same vector
      // stop the inner loop as soon as we hit identity so we avoid self interaction and duplicates
      if (&i == &j) {
        break;
      }

      auto dr = autopas::utils::ArrayMath::sub(i.getR(), j.getR());
      auto distanceSquare = autopas::utils::ArrayMath::dot(dr, dr);

      if (distanceSquare < _cutoff) {
        // first i then j because j has the smaller id
        _reference.emplace_back(j.getID(), i.getID());
      }
    }
  }
}

/**
 * Tests the Collision functor with every algorithm configuration in AutoPas.
 * The functor is applied, the colliding id's gathered and compared to the reference that was
 * computed in SetUpTestSuite().
 */
TEST_P(CollisionFunctorIntegrationTest, testAutoPasAlgorithm) {
  // This test does not really make sense if there are no close encounters in the dataset
  ASSERT_THAT(_reference, testing::Not(testing::IsEmpty()));

  const auto &[traversal, dataLayout, newton3, cellSizeFactor] = GetParam();

  // This is currently necessary until we implement the SoA functor
  if (dataLayout == autopas::DataLayoutOption::soa) {
    GTEST_SKIP_("SoAFunctor currently not implemented!");
  }

  CollisionFunctor functor(_cutoff, 10.0, 0.1 * _cutoff);

  // configure the AutoPas container
  autopas::AutoPas<Particle> autopas;
  // allow all container options since the traversal determines it uniquely
  autopas.setAllowedContainers({autopas::ContainerOption::getAllOptions()});
  autopas.setAllowedTraversals({traversal});
  autopas.setAllowedDataLayouts({dataLayout});
  autopas.setAllowedNewton3Options({newton3});
  autopas.setAllowedCellSizeFactors(autopas::NumberSetFinite{cellSizeFactor});
  autopas.setBoxMin(_boxMin);
  autopas.setBoxMax(_boxMax);
  autopas.setCutoff(_cutoff);
  autopas.init();

  // fill the container with the particles
  for (const auto &d : _debris) {
    autopas.addParticle(d);
  }

  try {
    autopas.iteratePairwise(&functor);
  } catch (const autopas::utils::ExceptionHandler::AutoPasException &e) {
    // There will be some tests generated that are invalid configurations but that is ok.
    GTEST_SKIP_(e.what());
  }

  // translate references to ids
  auto collisionPtrs = functor.getCollisions();
  std::vector<std::pair<size_t, size_t>> collisionIDs;
  collisionIDs.reserve(collisionPtrs.size());
  for (const auto &[pi, pj, dist] : collisionPtrs) {
    collisionIDs.emplace_back(pi->getID(), pj->getID());
  }

  EXPECT_THAT(collisionIDs, ::testing::UnorderedElementsAreArray(_reference));
}

// Generate tests for all configuration combinations
INSTANTIATE_TEST_SUITE_P(Generated,
                         CollisionFunctorIntegrationTest,
                         testing::Combine(testing::ValuesIn(autopas::TraversalOption::getAllOptions()),
                                          testing::ValuesIn(autopas::DataLayoutOption::getAllOptions()),
                                          testing::ValuesIn(autopas::Newton3Option::getAllOptions()),
                                          testing::ValuesIn({0.5, 1., 2.5})),
                         CollisionFunctorIntegrationTest::PrintToStringParamName());
