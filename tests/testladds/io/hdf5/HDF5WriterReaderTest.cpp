/**
 * @file HDF5WriterReaderTest.cpp
 * @author F. Gratl
 * @date 23.12.2021
 */

#include <gmock/gmock-matchers.h>

#include "HDF5WriterReaderTest.h"
#include "ladds/io/hdf5/HDF5Reader.h"
#include "ladds/io/hdf5/HDF5Writer.h"

#ifdef LADDS_HDF5
TEST_F(HDF5WriterReaderTest, WriteReadTestParticleData) {
  // 1. create some data
  constexpr size_t numParticles = 2;
  autopas::AutoPas<LADDS::Particle> autopas;
  autopas.setBoxMin({0., 0., 0.});
  autopas.setBoxMax(
      {static_cast<double>(numParticles), static_cast<double>(numParticles), static_cast<double>(numParticles)});
  autopas.init();

  std::vector<LADDS::Particle> particles;
  particles.reserve(numParticles);
  for (size_t i = 0; i < numParticles; ++i) {
    LADDS::Particle p{{static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)},
                      {1., 2., 3.},
                      i,
                      "dummy",
                      LADDS::Particle::ActivityState::evasive,
                      1.,
                      1.,
                      LADDS::Particle::calculateBcInv(0., 1., 1., 2.2)};
    autopas.addParticle(p);
    particles.push_back(p);
  }
  // sanity check
  ASSERT_GT(numParticles, 0);
  ASSERT_EQ(numParticles, autopas.getNumberOfParticles());
  ASSERT_EQ(numParticles, particles.size());

  // 2. write data
  constexpr auto filename = "WriteReadTestParticleData.h5";
  constexpr size_t iterationNr = 42;
  LADDS::HDF5Writer hdf5Writer(filename, true, 4);
  hdf5Writer.writeParticles(iterationNr, autopas);

  // 3. read data and check that read data is equal to generated data
  LADDS::HDF5Reader hdf5Reader(filename);
  EXPECT_THAT(hdf5Reader.readParticles(iterationNr, 2.2), ::testing::UnorderedElementsAreArray(particles))
      << "Particle data of initial iteration is not correct!";

  // 4. add more particles and check that they are added correctly to HDF5
  particles.emplace_back(std::array<double, 3>{0.5, 0.5, 0.5},
                         std::array<double, 3>{1., 2., 3.},
                         numParticles,
                         "dummy",
                         LADDS::Particle::ActivityState::evasivePreserving,
                         1.,
                         1.,
                         LADDS::Particle::calculateBcInv(0., 1., 1., 2.2),
                         std::numeric_limits<size_t>::max());
  autopas.addParticle(particles.back());

  // 5. write data
  hdf5Writer.writeParticles(iterationNr + 1, autopas);

  // 6. read data and check that read data is equal to generated data
  EXPECT_THAT(hdf5Reader.readParticles(iterationNr + 1, 2.2), ::testing::UnorderedElementsAreArray(particles))
      << "Particle data of second iteration is not correct!";

  // cleanup
  std::remove(filename);
}

TEST_F(HDF5WriterReaderTest, WriteReadTestCollisionData) {
  // 1. create some data
  constexpr size_t numParticles = 4;
  std::vector<LADDS::Particle> particles;
  particles.reserve(numParticles);
  for (size_t i = 0; i < numParticles; ++i) {
    particles.emplace_back<LADDS::Particle>({{0., 0., static_cast<double>(i)},
                                             {0., 0., 0.},
                                             i,
                                             "dummy",
                                             LADDS::Particle::ActivityState::passive,
                                             1.,
                                             1.,
                                             LADDS::Particle::calculateBcInv(0., 1., 1., 2.2)});
  }

  // These conjunctions are just randomly made up and have nothing to do with position data!
  LADDS::CollisionFunctor::CollisionCollectionT conjunctions;
  auto insertConjunction = [&](size_t idA, size_t idB) {
    const auto dr = autopas::utils::ArrayMath::sub(particles[idA].getR(), particles[idB].getR());
    const auto distanceSquare = autopas::utils::ArrayMath::dot(dr, dr);
    conjunctions.emplace_back(&particles[idA], &particles[idB], distanceSquare);
  };
  insertConjunction(1, 2);
  insertConjunction(2, 3);
  insertConjunction(2, 0);

  // 2. write data
  constexpr auto filename = "WriteReadTestCollisionData.h5";
  constexpr size_t iterationNr = 42;
  LADDS::HDF5Writer hdf5Writer(filename, true, 4);
  hdf5Writer.writeConjunctions(iterationNr, conjunctions);

  // 3. read data
  LADDS::HDF5Reader hdf5Reader(filename);
  auto conjunctionsHDF5 = hdf5Reader.readCollisions(iterationNr);

  // 4. check that read data is equal to generated data
  EXPECT_EQ(conjunctions.size(), conjunctionsHDF5.size());
  for (const auto &[ptrARef, ptrBRef, distRef] : conjunctions) {
    const auto idARef = static_cast<LADDS::HDF5Definitions::IntType>(ptrARef->getID());
    const auto idBRef = static_cast<LADDS::HDF5Definitions::IntType>(ptrBRef->getID());

    LADDS::HDF5Definitions::CollisionInfo collisionInfo{
        idARef, idBRef, static_cast<LADDS::HDF5Definitions::FloatType>(distRef)};
    EXPECT_THAT(conjunctionsHDF5, ::testing::Contains(collisionInfo));
  }
  // cleanup
  std::remove(filename);
}

/**
 * Write data to a HDF5 file, create a second writer that appends some more data, load and check that everything is
 * there.
 */
TEST_F(HDF5WriterReaderTest, AppendCheckpointTest) {
  // 1. create some data
  constexpr size_t numParticles = 4;
  AutoPas_t autopas;
  autopas.setBoxMin({0., 0., 0.});
  autopas.setBoxMax({10., 10., 10.});
  autopas.init();
  std::vector<LADDS::Particle> particles;
  particles.reserve(numParticles + 1);
  for (size_t i = 0; i < numParticles; ++i) {
    particles.emplace_back<LADDS::Particle>({{0., 0., static_cast<double>(i)},
                                             {0., 0., 0.},
                                             i,
                                             "dummy",
                                             LADDS::Particle::ActivityState::passive,
                                             1.,
                                             1.,
                                             LADDS::Particle::calculateBcInv(0., 1., 1., 2.2),
                                             std::numeric_limits<size_t>::max()});
    autopas.addParticle(particles.back());
  }

  // These conjunctions are just randomly made up and have nothing to do with position data!
  auto insertConjunction = [&](size_t idA, size_t idB, LADDS::CollisionFunctor::CollisionCollectionT &conjunctions) {
    const auto dr = autopas::utils::ArrayMath::sub(particles[idA].getR(), particles[idB].getR());
    const auto distanceSquare = autopas::utils::ArrayMath::dot(dr, dr);
    conjunctions.emplace_back(&particles[idA], &particles[idB], distanceSquare);
  };
  LADDS::CollisionFunctor::CollisionCollectionT conjunctionsStepA;
  insertConjunction(1, 2, conjunctionsStepA);
  insertConjunction(2, 3, conjunctionsStepA);

  // 2. write data (StepA)
  constexpr auto filename = "AppendCheckpointTest.h5";
  constexpr size_t iterationStepA{42};
  {
    LADDS::HDF5Writer hdf5WriterReplace(filename, true, 4);
    hdf5WriterReplace.writeParticles(iterationStepA, autopas);
    hdf5WriterReplace.writeConjunctions(iterationStepA, conjunctionsStepA);
  }

  // 3. new writer that appends (StepB)
  constexpr size_t iterationStepB{1337};
  LADDS::CollisionFunctor::CollisionCollectionT conjunctionsStepB;
  insertConjunction(2, 0, conjunctionsStepB);
  particles.emplace_back<LADDS::Particle>({{0., 0., static_cast<double>(numParticles)},
                                           {0., 0., 0.},
                                           numParticles,
                                           "dummy",
                                           LADDS::Particle::ActivityState::passive,
                                           1.,
                                           1.,
                                           LADDS::Particle::calculateBcInv(0., 1., 1., 2.2),
                                           std::numeric_limits<size_t>::max()});
  autopas.addParticle(particles.back());
  {
    LADDS::HDF5Writer hdf5WriterAppend(filename, false, 4);
    hdf5WriterAppend.writeParticles(iterationStepB, autopas);
    hdf5WriterAppend.writeConjunctions(iterationStepB, conjunctionsStepB);
  }

  // 4. check that all data is present
  LADDS::HDF5Reader hdf5Reader(filename);
  EXPECT_THAT(hdf5Reader.readParticles(iterationStepB, 2.2), ::testing::UnorderedElementsAreArray(particles))
      << "Particle data of of StepB is incorrect!";

  auto conjunctionsHDF5StepA = hdf5Reader.readCollisions(iterationStepA);
  EXPECT_EQ(conjunctionsStepA.size(), conjunctionsHDF5StepA.size());
  for (const auto &[ptrARef, ptrBRef, distRef] : conjunctionsStepA) {
    const auto idARef = static_cast<LADDS::HDF5Definitions::IntType>(ptrARef->getID());
    const auto idBRef = static_cast<LADDS::HDF5Definitions::IntType>(ptrBRef->getID());

    LADDS::HDF5Definitions::CollisionInfo collisionInfo{
        idARef, idBRef, static_cast<LADDS::HDF5Definitions::FloatType>(distRef)};
    EXPECT_THAT(conjunctionsHDF5StepA, ::testing::Contains(collisionInfo));
  }
  auto conjunctionsHDF5StepB = hdf5Reader.readCollisions(iterationStepB);
  EXPECT_EQ(conjunctionsStepB.size(), conjunctionsHDF5StepB.size());
  for (const auto &[ptrARef, ptrBRef, distRef] : conjunctionsStepB) {
    const auto idARef = static_cast<LADDS::HDF5Definitions::IntType>(ptrARef->getID());
    const auto idBRef = static_cast<LADDS::HDF5Definitions::IntType>(ptrBRef->getID());

    LADDS::HDF5Definitions::CollisionInfo collisionInfo{
        idARef, idBRef, static_cast<LADDS::HDF5Definitions::FloatType>(distRef)};
    EXPECT_THAT(conjunctionsHDF5StepB, ::testing::Contains(collisionInfo));
  }
}
#else
TEST(HDF5WriterReaderTest, TestWriterNotCompiledException) {
  EXPECT_THROW(HDF5Writer("foo", true, 0), std::runtime_error);
}
#endif
