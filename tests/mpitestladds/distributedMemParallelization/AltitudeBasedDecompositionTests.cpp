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
  config["sim"]["deltaT"] = 10.0;
  config["sim"]["maxAltitude"] = 85000.;
  config["sim"]["prop"]["coefficientOfDrag"] = 2.2;

  // optional parameters which are necessary for the tests here
  config["io"]["constellationCutoff"] = constellationCutoff;
  config["sim"]["collisionDistanceFactor"] = 1.;
  config["sim"]["iterations"] = 1;
  config["sim"]["minAltitude"] = 150.;
  config["sim"]["prop"]["useKEPComponent"] = false;
  config["sim"]["decompositionType"] = "Altitude";

  configReader = std::make_unique<LADDS::ConfigReader>(config, logger);

  decomposition = std::make_unique<LADDS::AltitudeBasedDecomposition>(*configReader);

  // Set ID ranges correctly in CFG to have unique IDs
  // (Can't do in constructor because we need the decomp to get the rank)
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(decomposition->getCommunicator(), &rank);
  const auto lengthIDRange = 1000;
  configReader->setNextSafeParticleID(lengthIDRange * rank);
  configReader->setLastSafeParticleID(lengthIDRange * (rank + 1) - 1);

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

  unsigned long numParticlesGlobal{};
  unsigned long numParticlesLocal = autopas->getNumberOfParticles();
  autopas::AutoPas_MPI_Allreduce(&numParticlesLocal,
                                 &numParticlesGlobal,
                                 1,
                                 AUTOPAS_MPI_UNSIGNED_LONG,
                                 AUTOPAS_MPI_SUM,
                                 decomposition->getCommunicator());

  ASSERT_EQ(numParticlesGlobal, 4ul) << "Expected 4 total particles but got " << numParticlesGlobal << " on rank "
                                     << rank;

  for (auto &particle : *autopas) {
    // Send particles from rank 0 to rank 1
    if (decomposition->getRank(particle.getPosition()) == 0) {
      std::cout << "Changing particle " << particle.getID() << " from rank " << rank << std::endl;
      particle.setPosition({0., 0., 7000});
    }
    // Send particles from rank 1 to rank 2
    else if (decomposition->getRank(particle.getPosition()) == 1) {
      std::cout << "Changing particle " << particle.getID() << " from rank " << rank << std::endl;
      particle.setPosition({0., 0., 7500});
    }
    // Send particles from rank 2 to rank 0
    else if (decomposition->getRank(particle.getPosition()) == 2) {
      std::cout << "Changing particle " << particle.getID() << " from rank " << rank << std::endl;
      particle.setPosition({0., 0., 6500});
    }
  }

  std::vector<LADDS::Particle> leavingParticles = decomposition->getAndRemoveLeavingParticles(*autopas);
  for (auto &particle : leavingParticles) {
    std::cout << "Rank " << rank << ": " << particle.getID() << " is leaving." << std::endl;
  }
  unsigned long numLeavingParticlesGlobal{};
  unsigned long numLeavingParticlesLocal = leavingParticles.size();
  autopas::AutoPas_MPI_Allreduce(&numLeavingParticlesLocal,
                                 &numLeavingParticlesGlobal,
                                 1,
                                 AUTOPAS_MPI_UNSIGNED_LONG,
                                 AUTOPAS_MPI_SUM,
                                 decomposition->getCommunicator());

  ASSERT_EQ(numLeavingParticlesGlobal, 3ul)
      << "Expected 3 leaving particles but got " << numLeavingParticlesGlobal << " on rank " << rank;

  std::array<size_t, 3> expectedIDs{0, 1000, 2000};
  for (auto &leavingParticle : leavingParticles) {
    std::cout << leavingParticle.getID() << std::endl;
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
 * Tests whether leaving particles in collisions are caught correctly
 * This test is expected to run on 8 MPI ranks
 */
TEST_F(AltitudeBasedDecompositionTests, testAltDecompCollisions) {
  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(decomposition->getCommunicator(), &numRanks);
  if (numRanks != 8) {
    GTEST_FAIL() << "Test is expected to be launched with 8 ranks but only has " << numRanks;
  }
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(decomposition->getCommunicator(), &rank);

  auto [csvWriter, accumulator, integrator] = simulation.initIntegrator(*autopas, *configReader);

  // initialize N particles on the local rank
  std::vector<LADDS::Particle> particles{};
  constexpr int N = 6;
  std::array<std::array<double, 3>, N> positions{{
      {0., 0., 6705.},
      {0., 0., 6705.5},
      {0., 0., 6707.},
      {0., 6705., 0},
      {0., 6707.5, 0},
      {0., 6707., 0},
  }};
  std::array<std::array<double, 3>, N> velocities{{
      {0., 0., 1.},
      {0., 0., 1.},
      {0., 0., 0.},
      {0., 0., 0.},
      {0., -1., 0.},
      {0., -1., 0.},
  }};
  for (size_t i = 0; i < N; ++i) {
    constexpr double mass = 42.;
    constexpr double radius = 42.;
    const LADDS::Particle p{
        positions[i],
        velocities[i],
        static_cast<size_t>(rank) * 1000ul + i,
        "Rank " + std::to_string(rank) + " p" + std::to_string(i),
        LADDS::Particle::ActivityState::passive,
        mass,
        radius,
        LADDS::Particle::calculateBcInv(std::numeric_limits<double>::quiet_NaN(), mass, radius, 2.2),
        std::numeric_limits<size_t>::max()};
    particles.push_back(p);
  }

  // std::cout << autopas::utils::ArrayUtils::to_string(decomposition->altitudeIntervals) << std::endl;

  // Add the particles we expect to collide
  LADDS::SatelliteLoader::addSatellitesToAutoPas(*autopas, particles, *decomposition, *configReader);

  for (auto &p : *autopas) {
    std::cout << "Rank " << rank << ": " << p.getID() << " at "
              << autopas::utils::ArrayUtils::to_string(p.getPosition()) << " was added." << std::endl;
  }

  // Check particles are in the right places
  if (rank < 2)
    ASSERT_EQ(autopas->getNumberOfParticles(), 3) << "Expected " << 3 << " particles on rank " << rank;
  else
    ASSERT_EQ(autopas->getNumberOfParticles(), 0) << "Expected " << 0 << " particles on rank " << rank;

  // Perform one timestep
  integrator->integrate(false);

  // (potentially) update the internal data structure and collect particles which are leaving the container.
  auto leavingParticles = autopas->updateContainer();

  // for altitude decomp we cannot rely on autopas square boxes, thus
  // we need to update the leaving particles manually as well
  auto manuallyLeavingParticles = decomposition->getAndRemoveLeavingParticles(*autopas);

  // add them to already leaving particles
  leavingParticles.insert(leavingParticles.end(), manuallyLeavingParticles.begin(), manuallyLeavingParticles.end());

  // check if the leaving particles are correct
  for (auto &p : leavingParticles)
    std::cout << "Rank " << rank << ": " << p.getID() << " at "
              << autopas::utils::ArrayUtils::to_string(p.getPosition()) << " is leaving." << std::endl;

  if (rank < 2)
    ASSERT_EQ(leavingParticles.size(), 2) << "Expected " << 2 << " particles on rank " << rank;
  else
    ASSERT_EQ(leavingParticles.size(), 0) << "Expected " << 0 << " particles on rank " << rank;

  auto [leaving_collisions, leaving_evasions] = LADDS::ParticleMigrationHandler::collisionDetectionAroundParticles(

      *autopas,
      leavingParticles,
      config["sim"]["deltaT"].as<double>() * autopas->getVerletRebuildFrequency(),
      8.,
      config["sim"]["collisionDistanceFactor"].as<double>(),
      0.05,
      0.1,
      true);

  if (rank < 2) std::cout << "Rank: " << rank << " Leaving collisions: " << leaving_collisions.size() << std::endl;

  // Check if the leaving particles collisions are correct
  if (rank < 2)
    ASSERT_EQ(leaving_collisions.size(), 1) << "Expected " << 0 << " particles on rank " << rank;
  else
    ASSERT_EQ(leaving_collisions.size(), 0) << "Expected " << 0 << " particles on rank " << rank;

  auto incomingParticles = decomposition->communicateParticles(leavingParticles, *autopas);

  for (auto &p : incomingParticles)
    std::cout << "Rank " << rank << ": " << p.getID() << " at "
              << autopas::utils::ArrayUtils::to_string(p.getPosition()) << " is arriving." << std::endl;

  // Check if the arriving particles are correct
  if (rank < 2)
    ASSERT_EQ(incomingParticles.size(), 2) << "Expected " << 2 << " particles on rank " << rank;
  else
    ASSERT_EQ(incomingParticles.size(), 0) << "Expected " << 0 << " particles on rank " << rank;

  auto [incoming_collisions, incoming_evasions] = LADDS::ParticleMigrationHandler::collisionDetectionAroundParticles(

      *autopas,
      incomingParticles,
      config["sim"]["deltaT"].as<double>() * autopas->getVerletRebuildFrequency(),
      8.,
      config["sim"]["collisionDistanceFactor"].as<double>(),
      0.1,
      0.05);

  if (rank < 2) std::cout << "Rank: " << rank << " Incoming Collisions: " << incoming_collisions.size() << std::endl;
  // Check if the incoming particles collisions are correct
  if (rank < 2)
    ASSERT_EQ(incoming_collisions.size(), 2) << "Expected " << 2 << " particles on rank " << rank;
  else
    ASSERT_EQ(incoming_collisions.size(), 0) << "Expected " << 0 << " particles on rank " << rank;

  // Print all collisions
  incoming_collisions.insert(incoming_collisions.end(), leaving_collisions.begin(), leaving_collisions.end());
  if (not incoming_collisions.empty()) {
    std::cout << "The following particles collided between ranks:" << std::endl;
    for (const auto &[p1, p2, _, __] : incoming_collisions) {
      std::cout << p1->getID() << " at " << autopas::utils::ArrayUtils::to_string(p1->getPosition()) << p2->getID()
                << " at " << autopas::utils::ArrayUtils::to_string(p2->getPosition()) << std::endl;
    }
  }
}