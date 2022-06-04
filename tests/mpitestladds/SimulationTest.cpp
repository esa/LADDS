/**
 * @file SimulationTest.cpp
 * @author F. Gratl
 * @date 02.06.22
 */

#include "SimulationTest.h"

/**
 * This test is expected to run on 8 MPI ranks
 */
TEST_F(SimulationTest, testRankMigration) {
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::mul;
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::sub;
  using autopas::utils::ArrayUtils::static_cast_array;

  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(decomposition->getCommunicator(), &numRanks);
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numRanks);
  if (numRanks != 8) {
    GTEST_FAIL() << "Test is expected to be launched with 8 ranks but only has " << numRanks;
  }
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(decomposition->getCommunicator(), &rank);

  configReader->setValue("sim/breakup/enabled", true);
  auto [csvWriter, accumulator, integrator] = simulation.initIntegrator(*autopas, *configReader);

  const auto localBoxMin = decomposition->getLocalBoxMin();
  const auto localBoxMax = decomposition->getLocalBoxMax();
  const auto localBoxLength = sub(localBoxMax, localBoxMin);
  const auto localBoxLengthHalf = mulScalar(localBoxLength, .5);
  const auto localBoxMid = add(localBoxMin, localBoxLengthHalf);

  // initialize N particles on the local rank
  std::vector<LADDS::Particle> particles;
  for (size_t i = 0; i < numRanks; ++i) {
    constexpr double mass = 42.;
    constexpr double radius = 42.;
    autopas->addParticle(
        {localBoxMid,
         {1., 2., 3.},
         static_cast<size_t>(rank) * 1000ul + i,
         "Rank " + std::to_string(rank) + " p" + std::to_string(i),
         LADDS::Particle::ActivityState::evasivePreserving,
         mass,
         radius,
         LADDS::Particle::calculateBcInv(std::numeric_limits<double>::quiet_NaN(), mass, radius, 2.2)});
  }
  std::cout << "RANK " << rank  << ": " << autopas::utils::ArrayUtils::to_string(localBoxMin) << " - " << autopas::utils::ArrayUtils::to_string(localBoxMax) << std::endl;
  // move the particles to other ranks
  auto particleIter = autopas->begin();
  for (int z = 0; z < 2; ++z) {
    for (int y = 0 ; y < 2; ++y) {
      for (int x = 0; x < 2; ++x) {
        const std::array<int, 3> rankGridIndex{x, y, z};
        const auto newPosition = add(localBoxLengthHalf, mul(localBoxLength, static_cast_array<double>(rankGridIndex)));
        // make sure all except one particle is moved away
        ASSERT_TRUE(
            (x == 0 and y == 0 and z == 0) or
            autopas::utils::notInBox(newPosition, decomposition->getLocalBoxMin(), decomposition->getLocalBoxMax()))
            << "Updated positions should not be in the same box anymore!";
        particleIter->setPosition(newPosition);
        std::cout << "RANK " << rank  << ": " << particleIter->getIdentifier() << " newPos: " << autopas::utils::ArrayUtils::to_string(newPosition) << std::endl;
      }
    }
  }
  // migrate particles
  auto leavingParticles = autopas->updateContainer();
  ASSERT_EQ(autopas->getNumberOfParticles(), 1) << "Expected all except one particle to have emigrated";
  simulation.communicateParticles(leavingParticles, *autopas, *decomposition);
  EXPECT_EQ(autopas->getNumberOfParticles(), numRanks) << "Not enough particles immigrated! " << [&]() {
    std::vector<std::string> localParticleIdentifiers;
    for (const auto &p : *autopas) {
      localParticleIdentifiers.push_back(p.getIdentifier());
    }
    return autopas::utils::ArrayUtils::to_string(localParticleIdentifiers);
  }();
}