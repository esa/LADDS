/**
 * @file RegularGridDecompositionTests.cpp
 * @author F. Gratl
 * @date 02.06.22
 */

#include "RegularGridDecompositionTests.h"

#include "autopas/utils/WrapOpenMP.h"
#include "ladds/TypeDefinitions.h"
#include "ladds/distributedMemParallelization/ParticleMigrationHandler.h"
#include "ladds/io/SatelliteLoader.h"

RegularGridDecompositionTests::RegularGridDecompositionTests()
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
  config["sim"]["decompositionType"] = "RegularGrid";

  configReader = std::make_unique<LADDS::ConfigReader>(config, logger);

  decomposition = std::make_unique<LADDS::RegularGridDecomposition>(*configReader);

  autopas = simulation.initAutoPas(*configReader, *decomposition);
}

RegularGridDecompositionTests::~RegularGridDecompositionTests() {
  // reset omp max threads
  autopas::autopas_set_num_threads(maxThreadsBefore);
}

/**
 * This test is expected to run on 8 MPI ranks
 */
TEST_F(RegularGridDecompositionTests, testParticleMigrationHandler) {
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
  const auto incomingParticles = decomposition->communicateParticles(leavingParticles, *autopas);
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

/**
 * Tests whether leaving particles in collisions are caught correctly
 * This test is expected to run on 8 MPI ranks
 */
TEST_F(RegularGridDecompositionTests, testGridDecompCollisions) {
  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(decomposition->getCommunicator(), &numRanks);
  if (numRanks != 8) {
    GTEST_FAIL() << "Test is expected to be launched with 8 ranks but only has " << numRanks;
  }
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(decomposition->getCommunicator(), &rank);

  // otherwise we create A LOT of acceleration in this test :)
  config["sim"]["prop"]["useKEPComponent"] = false;

  auto [csvWriter, accumulator, integrator] = simulation.initIntegrator(*autopas, *configReader);

  // initialize N particles on the local rank
  std::vector<LADDS::Particle> particles{};
  constexpr int N = 3;
  std::array<std::array<double, 3>, N> positions{{
      {0., 0., -1.},
      {0., 0., -0.5},
      {0., 0., 4.},
  }};
  std::array<std::array<double, 3>, N> velocities{{
      {0., 0., 5.},
      {0., 0., 5.},
      {0., 0., 0.},
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

  // auto localMin = decomposition->getLocalBoxMin();
  // auto localMax = decomposition->getLocalBoxMax();
  // std::cout << "Rank " << rank << " - min: " << autopas::utils::ArrayUtils::to_string(localMin)
  //           << " max: " << autopas::utils::ArrayUtils::to_string(localMax) << std::endl;

  // Add the particles we expect to collide
  LADDS::SatelliteLoader::addSatellitesToAutoPas(*autopas, particles, *decomposition, *configReader);

  for (auto &p : *autopas) {
    std::cout << "Rank " << rank << ": " << p.getID() << " at "
              << autopas::utils::ArrayUtils::to_string(p.getPosition()) << " was added." << std::endl;
  }

  // Check particles are in the right places
  if (rank == 6)
    ASSERT_EQ(autopas->getNumberOfParticles(), 2) << "Expected " << 3 << " particles on rank " << rank;
  else if (rank == 7)
    ASSERT_EQ(autopas->getNumberOfParticles(), 1) << "Expected " << 2 << " particles on rank " << rank;
  else
    ASSERT_EQ(autopas->getNumberOfParticles(), 0) << "Expected " << 0 << " particles on rank " << rank;

  // Perform one timestep
  integrator->integrate(false);

  // (potentially) update the internal data structure and collect particles which are leaving the container.
  auto leavingParticles = autopas->updateContainer();

  // check if the leaving particles are correct
  for (auto &p : leavingParticles)
    std::cout << "Rank " << rank << ": " << p.getID() << " at "
              << autopas::utils::ArrayUtils::to_string(p.getPosition()) << " is leaving." << std::endl;

  if (rank == 6)
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

  if (rank == 6 or rank == 7)
    std::cout << "Rank: " << rank << " Leaving collisions: " << leaving_collisions.size() << std::endl;

  // Check if the leaving particles collisions are correct
  if (rank == 6)
    ASSERT_EQ(leaving_collisions.size(), 1) << "Expected " << 1 << " particles on rank " << rank;
  else
    ASSERT_EQ(leaving_collisions.size(), 0) << "Expected " << 0 << " particles on rank " << rank;

  auto incomingParticles = decomposition->communicateParticles(leavingParticles, *autopas);

  for (auto &p : incomingParticles)
    std::cout << "Rank " << rank << ": " << p.getID() << " at "
              << autopas::utils::ArrayUtils::to_string(p.getPosition()) << " is arriving." << std::endl;

  // Check if the arriving particles are correct
  // https://github.com/esa/LADDS/issues/152
  GTEST_SKIP_("RegularGridDecomposition currently not operational!");
  if (rank == 7)
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
      0.05,
      false);

  if (rank == 6 or rank == 7)
    std::cout << "Rank: " << rank << " Incoming Collisions: " << incoming_collisions.size() << std::endl;
  // Check if the incoming particles collisions are correct
  if (rank == 6)
    ASSERT_EQ(incoming_collisions.size(), 0) << "Expected " << 0 << " particles on rank " << rank;
  else
    ASSERT_EQ(incoming_collisions.size(), 0) << "Expected " << 0 << " particles on rank " << rank;

  // Print all collisions
  incoming_collisions.insert(incoming_collisions.end(), leaving_collisions.begin(), leaving_collisions.end());
  if (not incoming_collisions.empty()) {
    std::cout << "The following particles collided between ranks:" << std::endl;
    for (const auto &[p1, p2, _, __] : incoming_collisions) {
      std::cout << p1->getID() << " at " << autopas::utils::ArrayUtils::to_string(p1->getPosition()) << " "
                << p2->getID() << " at " << autopas::utils::ArrayUtils::to_string(p2->getPosition()) << std::endl;
    }
  }
}