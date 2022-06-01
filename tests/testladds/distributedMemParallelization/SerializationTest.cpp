/**
 * @file SerializationTest.cpp
 * @author F. Gratl
 * @date 31/05/2022
 */

#include "SerializationTest.h"

#include "ladds/distributedMemParallelization/Serialization.h"
#include "ladds/particle/Particle.h"

TEST(SerializationTest, serializeAndDeserializeTest) {
  std::vector<LADDS::Particle> originalParticles{{{1., 2., 3.},
                                                  {4., 5., 6.},
                                                  42ul,
                                                  "2022-012XYZ",
                                                  LADDS::Particle::ActivityState::evasivePreserving,
                                                  4.2,
                                                  13.37,
                                                  0.5},
                                                 {{10., 20., 30.},
                                                  {40., 50., 60.},
                                                  420ul,
                                                  "3022-012XYZ",
                                                  LADDS::Particle::ActivityState::evasive,
                                                  40.2,
                                                  130.37,
                                                  0.05}};
  originalParticles[0].setAccT0({7., 8., 9.});
  originalParticles[0].setAccT1({10., 11., 12.});
  originalParticles[1].setAccT0({70., 80., 90.});
  originalParticles[1].setAccT1({100., 110., 120.});

  std::vector<char> serializedParticle;
  LADDS::Serialization::serializeParticles(originalParticles.begin(), originalParticles.end(), serializedParticle);
  std::vector<LADDS::Particle> deserializedParticles{};
  LADDS::Serialization::deserializeParticles(serializedParticle, deserializedParticles);

  for (int i = 0; i < originalParticles.size(); ++i) {
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::posX>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::posX>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::posY>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::posY>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::posZ>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::posZ>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::forceX>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::forceX>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::forceY>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::forceY>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::forceZ>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::forceZ>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::velocityX>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::velocityX>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::velocityY>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::velocityY>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::velocityZ>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::velocityZ>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::ownershipState>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::ownershipState>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::acc_t0X>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::acc_t0X>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::acc_t0Y>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::acc_t0Y>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::acc_t0Z>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::acc_t0Z>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::acc_t1X>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::acc_t1X>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::acc_t1Y>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::acc_t1Y>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::acc_t1Z>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::acc_t1Z>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::aom>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::aom>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::mass>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::mass>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::radius>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::radius>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::bc_inv>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::bc_inv>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::activityState>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::activityState>())
        << "for i=" << i;
    EXPECT_EQ(deserializedParticles[i].get<LADDS::Particle::AttributeNames::identifier>(),
              originalParticles[i].get<LADDS::Particle::AttributeNames::identifier>())
        << "for i=" << i;
  }
}