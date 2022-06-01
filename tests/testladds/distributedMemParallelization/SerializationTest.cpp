/**
 * @file SerializationTest.cpp
 * @author F. Gratl
 * @date 31/05/2022
 */

#include "SerializationTest.h"

#include "ladds/distributedMemParallelization/Serialization.h"
#include "ladds/particle/Particle.h"

TEST(SerializationTest, serializeAndDeserializeTest) {
  LADDS::Particle originalParticle({1., 2., 3.},
                                   {4., 5., 6.},
                                   42ul,
                                   "2022-012XYZ",
                                   LADDS::Particle::ActivityState::evasivePreserving,
                                   4.2,
                                   13.37,
                                   0.5);
  originalParticle.setAccT0({7., 8., 9.});
  originalParticle.setAccT1({10., 11., 12.});

  std::vector<char> serializedParticle;
  LADDS::Serialization::serializeParticle(originalParticle, serializedParticle);
  std::vector<LADDS::Particle> deserializedParticles{};
  LADDS::Serialization::deserializeParticles(serializedParticle, deserializedParticles);

  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::posX>(),
            originalParticle.get<LADDS::Particle::AttributeNames::posX>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::posY>(),
            originalParticle.get<LADDS::Particle::AttributeNames::posY>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::posZ>(),
            originalParticle.get<LADDS::Particle::AttributeNames::posZ>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::forceX>(),
            originalParticle.get<LADDS::Particle::AttributeNames::forceX>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::forceY>(),
            originalParticle.get<LADDS::Particle::AttributeNames::forceY>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::forceZ>(),
            originalParticle.get<LADDS::Particle::AttributeNames::forceZ>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::velocityX>(),
            originalParticle.get<LADDS::Particle::AttributeNames::velocityX>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::velocityY>(),
            originalParticle.get<LADDS::Particle::AttributeNames::velocityY>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::velocityZ>(),
            originalParticle.get<LADDS::Particle::AttributeNames::velocityZ>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::ownershipState>(),
            originalParticle.get<LADDS::Particle::AttributeNames::ownershipState>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::acc_t0X>(),
            originalParticle.get<LADDS::Particle::AttributeNames::acc_t0X>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::acc_t0Y>(),
            originalParticle.get<LADDS::Particle::AttributeNames::acc_t0Y>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::acc_t0Z>(),
            originalParticle.get<LADDS::Particle::AttributeNames::acc_t0Z>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::acc_t1X>(),
            originalParticle.get<LADDS::Particle::AttributeNames::acc_t1X>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::acc_t1Y>(),
            originalParticle.get<LADDS::Particle::AttributeNames::acc_t1Y>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::acc_t1Z>(),
            originalParticle.get<LADDS::Particle::AttributeNames::acc_t1Z>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::aom>(),
            originalParticle.get<LADDS::Particle::AttributeNames::aom>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::mass>(),
            originalParticle.get<LADDS::Particle::AttributeNames::mass>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::radius>(),
            originalParticle.get<LADDS::Particle::AttributeNames::radius>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::bc_inv>(),
            originalParticle.get<LADDS::Particle::AttributeNames::bc_inv>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::activityState>(),
            originalParticle.get<LADDS::Particle::AttributeNames::activityState>());
  EXPECT_EQ(deserializedParticles[0].get<LADDS::Particle::AttributeNames::identifier>(),
            originalParticle.get<LADDS::Particle::AttributeNames::identifier>());
}