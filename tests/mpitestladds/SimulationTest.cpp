/**
 * @file SimulationTest.cpp
 * @author F. Gratl
 * @date 02.06.22
 */

#include "SimulationTest.h"

#include "autopas/utils/WrapOpenMP.h"

SimulationTest::SimulationTest()
    : maxThreadsBefore(autopas::autopas_get_max_threads()), logger("SimulationTestLogger"), simulation(logger) {
  // make sure to only use one thread
  autopas::autopas_set_num_threads(1);

  logger.get()->set_level(LADDS::Logger::Level::err);

  // initialize a minimal default configuration
  config["autopas"]["cutoff"] = 80.;
  config["sim"]["breakup"]["enabled"] = false;
  config["sim"]["deltaT"] = 1.0;
  config["sim"]["maxAltitude"] = 85000.;
  config["sim"]["prop"]["coefficientOfDrag"] = 2.2;

  // optional parameters which are necessary for the tests here
  config["io"]["constellationCutoff"] = constellationCutoff;
  config["sim"]["collisionDistanceFactor"] = 1.;
  config["sim"]["iterations"] = 1;
  config["sim"]["minAltitude"] = 150.;
  config["sim"]["prop"]["useKEPComponent"] = true;

  configReader = std::make_unique<LADDS::ConfigReader>(config, logger);

  decomposition = std::make_unique<LADDS::RegularGridDecomposition>(*configReader);

  autopas = simulation.initAutoPas(*configReader, *decomposition);
}

SimulationTest::~SimulationTest() {
  // reset omp max threads
  autopas::autopas_set_num_threads(maxThreadsBefore);
}

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
  if (numRanks != 8) {
    GTEST_FAIL() << "Test is expected to be launched with 8 ranks but only has " << numRanks;
  }
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(decomposition->getCommunicator(), &rank);

  configReader->setValue("sim/breakup/enabled", true);
  auto [csvWriter, accumulator, integrator] = simulation.initIntegrator(*autopas, *configReader);

  const auto globalBoxMin = decomposition->getGlobalBoxMin();
  const auto localBoxMin = decomposition->getLocalBoxMin();
  const auto localBoxMax = decomposition->getLocalBoxMax();
  const auto localBoxLength = sub(localBoxMax, localBoxMin);
  const auto localBoxLengthHalf = mulScalar(localBoxLength, .5);
  const auto localBoxMid = add(localBoxMin, localBoxLengthHalf);

  // initialize N particles on the local rank
  std::vector<LADDS::Particle> particles{};
  for (size_t i = 0; i < numRanks; ++i) {
    constexpr double mass = 42.;
    constexpr double radius = 42.;
    const LADDS::Particle p{
        localBoxMid,
        {1., 2., 3.},
        static_cast<size_t>(rank) * 1000ul + i,
        "Rank " + std::to_string(rank) + " p" + std::to_string(i),
        LADDS::Particle::ActivityState::evasivePreserving,
        mass,
        radius,
        LADDS::Particle::calculateBcInv(std::numeric_limits<double>::quiet_NaN(), mass, radius, 2.2),
        std::numeric_limits<size_t>::max()};
    autopas->addParticle(p);
  }
  // move the particles to other ranks
  auto particleIter = autopas->begin();
  for (int x = 0; x < 2; ++x) {
    for (int y = 0; y < 2; ++y) {
      for (int z = 0; z < 2; ++z) {
        const std::array<int, 3> rankGridIndex{x, y, z};
        // newPos = globalBoxMin + (localBoxLength / 2) + (localBoxLength * index)
        const auto newPosition =
            add(globalBoxMin, add(localBoxLengthHalf, mul(localBoxLength, static_cast_array<double>(rankGridIndex))));
        // make sure all except one particle is moved away
        ASSERT_TRUE(
            (x == rankGridIndex[0] and y == rankGridIndex[1] and z == rankGridIndex[2]) or
            autopas::utils::notInBox(newPosition, decomposition->getLocalBoxMin(), decomposition->getLocalBoxMax()))
            << "Updated positions should not be in the same box anymore!\n"
            << "Box: " << autopas::utils::ArrayUtils::to_string(decomposition->getLocalBoxMin()) << " - "
            << autopas::utils::ArrayUtils::to_string(decomposition->getLocalBoxMax()) << "\n"
            << particleIter->toString();
        particleIter->setPosition(newPosition);
        ++particleIter;
      }
    }
  }
  // migrate particles
  auto leavingParticles = autopas->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 7) << "Expected all except one particle to have left.";
  ASSERT_EQ(autopas->getNumberOfParticles(), 1) << "Expected exactly one particle to remain.";
  const auto incomingParticles = simulation.communicateParticles(leavingParticles, *autopas, *decomposition);
  EXPECT_EQ(incomingParticles.size(), numRanks - 1);
  for (const auto &p : incomingParticles) {
    autopas->addParticle(p);
  }
  EXPECT_EQ(autopas->getNumberOfParticles(), numRanks)
      << "Rank " << rank << ": Wrong number of particles immigrated! " << [&]() {
           std::vector<std::string> localParticleIdentifiers;
           for (const auto &p : *autopas) {
             localParticleIdentifiers.push_back(p.getIdentifier());
           }
           return autopas::utils::ArrayUtils::to_string(localParticleIdentifiers);
         }();
}