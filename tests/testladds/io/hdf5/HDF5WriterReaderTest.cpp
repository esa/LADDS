/**
 * @file HDF5WriterReaderTest.cpp
 * @author F. Gratl
 * @date 23.12.2021
 */

#include "HDF5WriterReaderTest.h"

#include <gmock/gmock-matchers.h>

#include "ladds/io/hdf5/HDF5Reader.h"
#include "ladds/io/hdf5/HDF5Writer.h"

#ifdef LADDS_HDF5
TEST_F(HDF5WriterReaderTest, WriteReadTestParticleData) {
  // 1. create some data
  constexpr size_t numParticles = 10;
  autopas::AutoPas<Particle> autopas;
  autopas.setBoxMin({0., 0., 0.});
  autopas.setBoxMax(
      {static_cast<double>(numParticles), static_cast<double>(numParticles), static_cast<double>(numParticles)});
  autopas.init();

  std::vector<Particle> particles;
  particles.reserve(numParticles);
  for (size_t i = 0; i < numParticles; ++i) {
    Particle p{{static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)},
               {1., 2., 3.},
               i,
               Particle::ActivityState::evasive,
               1.,
               1.};
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
  HDF5Writer hdf5Writer(filename, 4);
  hdf5Writer.writeParticles(iterationNr, autopas);

  // 3. read data and check that read data is equal to generated data
  HDF5Reader hdf5Reader(filename);
  EXPECT_THAT(hdf5Reader.readParticles(iterationNr), ::testing::UnorderedElementsAreArray(particles))
      << "Particle data of initial iteration is not correct!";

  // 4. add more particles and check that they are added correctly to HDF5
  particles.emplace_back(std::array<double, 3>{0.5, 0.5, 0.5},
                         std::array<double, 3>{1., 2., 3.},
                         numParticles,
                         Particle::ActivityState::evasivePreserving,
                         1.,
                         1.);
  autopas.addParticle(particles.back());

  // 5. write data
  hdf5Writer.writeParticles(iterationNr + 1, autopas);

  // 6. read data and check that read data is equal to generated data
  EXPECT_THAT(hdf5Reader.readParticles(iterationNr + 1), ::testing::UnorderedElementsAreArray(particles))
      << "Particle data of second iteration is not correct!";
  ;

  // cleanup
  std::remove(filename);
}

TEST_F(HDF5WriterReaderTest, WriteReadTestCollisionData) {
  // 1. create some data
  constexpr size_t numParticles = 4;
  std::vector<Particle> particles;
  particles.reserve(numParticles);
  for (size_t i = 0; i < numParticles; ++i) {
    particles.emplace_back<Particle>(
        {{0., 0., static_cast<double>(i)}, {0., 0., 0.}, i, Particle::ActivityState::passive, 1., 1.});
  }

  // These conjunctions are just randomly made up and have nothing to do with position data!
  CollisionFunctor::CollisionCollectionT conjunctions;
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
  HDF5Writer hdf5Writer(filename, 4);
  hdf5Writer.writeConjunctions(iterationNr, conjunctions);

  // 3. read data
  HDF5Reader hdf5Reader(filename);
  auto conjunctionsHDF5 = hdf5Reader.readCollisions(iterationNr);

  // 4. check that read data is equal to generated data
  EXPECT_EQ(conjunctions.size(), conjunctionsHDF5.size());
  for (const auto &[ptrARef, ptrBRef, distRef] : conjunctions) {
    const auto idARef = static_cast<HDF5Definitions::IntType>(ptrARef->getID());
    const auto idBRef = static_cast<HDF5Definitions::IntType>(ptrBRef->getID());

    HDF5Definitions::CollisionInfo collisionInfo{idARef, idBRef, static_cast<HDF5Definitions::FloatType>(distRef)};
    EXPECT_THAT(conjunctionsHDF5, ::testing::Contains(collisionInfo));
  }
  // cleanup
  std::remove(filename);
}
#else
TEST(HDF5WriterReaderTest, TestReaderNotCompiledException) {
  EXPECT_THROW(HDF5Reader("foo"), std::runtime_error);
}
TEST(HDF5WriterReaderTest, TestWriterNotCompiledException) {
  EXPECT_THROW(HDF5Writer("foo", 0), std::runtime_error);
}
#endif