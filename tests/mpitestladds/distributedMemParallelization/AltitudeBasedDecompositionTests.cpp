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

  AutoPas_MPI_Barrier(decomposition->getCommunicator());

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