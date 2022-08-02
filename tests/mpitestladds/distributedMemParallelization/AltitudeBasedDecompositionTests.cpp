/**
 * @file AltitudeBasedDecompositionTests.cpp
 * @author P. Gomez
 * @date 2022-08-01
 */

#include "AltitudeBasedDecompositionTests.h"

#include "autopas/utils/WrapOpenMP.h"
#include "ladds/TypeDefinitions.h"
#include "ladds/distributedMemParallelization/ParticleMigrationHandler.h"
#include "ladds/io/SatelliteLoader.h"
#include "satellitePropagator/physics/Constants.h"

AltitudeBasedDecompositionTests::AltitudeBasedDecompositionTests()
    : maxThreadsBefore(autopas::autopas_get_max_threads()), logger(LADDS_SPD_LOGGER_NAME), simulation(logger) {
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
  config["sim"]["decompositionType"] = "Altitude";

  configReader = std::make_unique<LADDS::ConfigReader>(config, logger);

  decomposition = std::make_unique<LADDS::AltitudeBasedDecomposition>(*configReader);

  autopas = simulation.initAutoPas(*configReader, *decomposition);
}

AltitudeBasedDecompositionTests::~AltitudeBasedDecompositionTests() {
  // reset omp max threads
  autopas::autopas_set_num_threads(maxThreadsBefore);
}

/**
 * Tests whether the decomposition of altitudes is correct.
 * This test is expected to run on 8 MPI ranks
 */
TEST_F(AltitudeBasedDecompositionTests, testAltitudeIntervalBalancing) {
  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(decomposition->getCommunicator(), &numRanks);
  if (numRanks != 8) {
    GTEST_FAIL() << "Test is expected to be launched with 8 ranks but only has " << numRanks;
  }
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(decomposition->getCommunicator(), &rank);

  // initialize N particles on the local rank
  std::vector<LADDS::Particle> particles{};
  // Number of particles
  int N = 64;
  for (size_t i = 0; i < N; ++i) {
    constexpr double mass = 42.;
    constexpr double radius = 42.;
    const double altitude = Physics::R_EARTH + 175 + 10 * i;
    const LADDS::Particle p{
        {0., 0., altitude},
        {1., 2., 3.},
        static_cast<size_t>(rank) * 1000ul + i,
        "Rank " + std::to_string(rank) + " p" + std::to_string(i),
        LADDS::Particle::ActivityState::evasivePreserving,
        mass,
        radius,
        LADDS::Particle::calculateBcInv(std::numeric_limits<double>::quiet_NaN(), mass, radius, 2.2),
        std::numeric_limits<size_t>::max()};
    particles.push_back(p);
  }

  LADDS::SatelliteLoader::addSatellitesToAutoPas(*autopas, particles, *decomposition, *configReader);

  ASSERT_TRUE(autopas->getNumberOfParticles() <= static_cast<int>((N / 8.)) + 1)
      << "Expected roughly equal distribution of particles per rank.";
  ASSERT_TRUE(autopas->getNumberOfParticles() >= static_cast<int>((N / 8.)) - 1)
      << "Expected roughly equal distribution of particles per rank.";
}

/**
 * Tests whether leaving particles are identified correctly
 * This test is expected to run on 8 MPI ranks
 */
TEST_F(AltitudeBasedDecompositionTests, testAltDecompParticleMigration) {
  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(decomposition->getCommunicator(), &numRanks);
  if (numRanks != 8) {
    GTEST_FAIL() << "Test is expected to be launched with 8 ranks but only has " << numRanks;
  }
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(decomposition->getCommunicator(), &rank);

  // initialize N particles on the local rank
  std::vector<LADDS::Particle> particles{};
  // Number of particles, this time less than ranks
  int N = 4;
  for (size_t i = 0; i < N; ++i) {
    constexpr double mass = 42.;
    constexpr double radius = 42.;
    const double altitude = Physics::R_EARTH + 175 + 500 * i;
    const LADDS::Particle p{
        {0., 0., altitude},
        {1., 2., 3.},
        static_cast<size_t>(rank) * 1000ul + i,
        "Rank " + std::to_string(rank) + " p" + std::to_string(i),
        LADDS::Particle::ActivityState::evasivePreserving,
        mass,
        radius,
        LADDS::Particle::calculateBcInv(std::numeric_limits<double>::quiet_NaN(), mass, radius, 2.2),
        std::numeric_limits<size_t>::max()};
    particles.push_back(p);
  }

  LADDS::SatelliteLoader::addSatellitesToAutoPas(*autopas, particles, *decomposition, *configReader);

  std::cout << autopas::utils::ArrayUtils::to_string(decomposition->altitudeIntervals) << std::endl;

  for (auto &particle : *autopas) {
    std::cout << particle << std::endl;
    // Send particles from rank 0 to rank 1
    if (decomposition->getRank(particle.getPosition()) == 0) particle.setPosition({0., 0., 7000});
    // Send particles from rank 1 to rank 2
    if (decomposition->getRank(particle.getPosition()) == 1) particle.setPosition({0., 0., 7500});
    // Send particles from rank 2 to rank 0
    if (decomposition->getRank(particle.getPosition()) == 2) particle.setPosition({0., 0., 6500});
  }

  std::vector<LADDS::Particle> leavingParticles = decomposition->getAndRemoveLeavingParticles(*autopas);
  ASSERT_EQ(leavingParticles.size(), 3ul) << "Expected 3 leaving particles.";

  std::array<size_t, 3> expectedIDs{0, 1, 2};
  for (auto &leavingParticle : leavingParticles) {
    ASSERT_TRUE(std::find(expectedIDs.begin(), expectedIDs.end(), leavingParticle.getID()) != expectedIDs.end())
        << "Expected leaving particle with ID " << leavingParticle.getID() << " to be in expectedIDs.";
  }

  decomposition->communicateParticles(leavingParticles, *autopas);

  for (auto &p : leavingParticles) {
    if (p.getID() == 0) ASSERT_TRUE(rank == 1) << "Expected leaving particles to be on rank 1.";
    if (p.getID() == 1) ASSERT_TRUE(rank == 2) << "Expected leaving particles to be on rank 2.";
    if (p.getID() == 2) ASSERT_TRUE(rank == 0) << "Expected leaving particles to be on rank 0.";
  }
}

/**
 * This test is expected to run on 8 MPI ranks
 */
TEST_F(AltitudeBasedDecompositionTests, testParticleMigration) {
  // using autopas::utils::ArrayMath::add;
  // using autopas::utils::ArrayMath::mul;
  // using autopas::utils::ArrayMath::mulScalar;
  // using autopas::utils::ArrayMath::sub;
  // using autopas::utils::ArrayUtils::static_cast_array;

  // int numRanks{};
  // autopas::AutoPas_MPI_Comm_size(decomposition->getCommunicator(), &numRanks);
  // if (numRanks != 8) {
  //   GTEST_FAIL() << "Test is expected to be launched with 8 ranks but only has " << numRanks;
  // }
  // int rank{};
  // autopas::AutoPas_MPI_Comm_rank(decomposition->getCommunicator(), &rank);

  // configReader->setValue("sim/breakup/enabled", true);
  // auto [csvWriter, accumulator, integrator] = simulation.initIntegrator(*autopas, *configReader);

  // const auto globalBoxMin = decomposition->getGlobalBoxMin();
  // const auto localBoxMin = decomposition->getLocalBoxMin();
  // const auto localBoxMax = decomposition->getLocalBoxMax();
  // const auto localBoxLength = sub(localBoxMax, localBoxMin);
  // const auto localBoxLengthHalf = mulScalar(localBoxLength, .5);
  // const auto localBoxMid = add(localBoxMin, localBoxLengthHalf);

  // // initialize N particles on the local rank
  // std::vector<LADDS::Particle> particles{};
  // for (size_t i = 0; i < numRanks; ++i) {
  //   constexpr double mass = 42.;
  //   constexpr double radius = 42.;
  //   const LADDS::Particle p{
  //       localBoxMid,
  //       {1., 2., 3.},
  //       static_cast<size_t>(rank) * 1000ul + i,
  //       "Rank " + std::to_string(rank) + " p" + std::to_string(i),
  //       LADDS::Particle::ActivityState::evasivePreserving,
  //       mass,
  //       radius,
  //       LADDS::Particle::calculateBcInv(std::numeric_limits<double>::quiet_NaN(), mass, radius, 2.2),
  //       std::numeric_limits<size_t>::max()};
  //   autopas->addParticle(p);
  // }
  // // move the particles to other ranks
  // auto particleIter = autopas->begin();
  // for (int x = 0; x < 2; ++x) {
  //   for (int y = 0; y < 2; ++y) {
  //     for (int z = 0; z < 2; ++z) {
  //       const std::array<int, 3> rankGridIndex{x, y, z};
  //       // newPos = globalBoxMin + (localBoxLength / 2) + (localBoxLength * index)
  //       const auto newPosition =
  //           add(globalBoxMin, add(localBoxLengthHalf, mul(localBoxLength,
  //           static_cast_array<double>(rankGridIndex))));
  //       // make sure all except one particle is moved away
  //       ASSERT_TRUE(
  //           (x == rankGridIndex[0] and y == rankGridIndex[1] and z == rankGridIndex[2]) or
  //           autopas::utils::notInBox(newPosition, decomposition->getLocalBoxMin(), decomposition->getLocalBoxMax()))
  //           << "Updated positions should not be in the same box anymore!\n"
  //           << "Box: " << autopas::utils::ArrayUtils::to_string(decomposition->getLocalBoxMin()) << " - "
  //           << autopas::utils::ArrayUtils::to_string(decomposition->getLocalBoxMax()) << "\n"
  //           << particleIter->toString();
  //       particleIter->setPosition(newPosition);
  //       ++particleIter;
  //     }
  //   }
  // }
  // // migrate particles
  // auto leavingParticles = autopas->updateContainer();
  // EXPECT_EQ(leavingParticles.size(), 7) << "Expected all except one particle to have left.";
  // ASSERT_EQ(autopas->getNumberOfParticles(), 1) << "Expected exactly one particle to remain.";
  // const auto incomingParticles = decomposition->communicateParticles(leavingParticles, *autopas);
  // EXPECT_EQ(incomingParticles.size(), numRanks - 1);
  // for (const auto &p : incomingParticles) {
  //   autopas->addParticle(p);
  // }
  // EXPECT_EQ(autopas->getNumberOfParticles(), numRanks)
  //     << "Rank " << rank << ": Wrong number of particles immigrated! " << [&]() {
  //          std::vector<std::string> localParticleIdentifiers;
  //          for (const auto &p : *autopas) {
  //            localParticleIdentifiers.push_back(p.getIdentifier());
  //          }
  //          return autopas::utils::ArrayUtils::to_string(localParticleIdentifiers);
  //        }();
}