#include "SatelliteToParticleConverterTest.h"

#include <breakupModel/model/Satellite.h>
#include <gmock/gmock-matchers.h>

#include "ladds/particle/SatelliteToParticleConverter.h"
/**
 * Sets up a particle and tries to convert it to a satellites.
 *
 */
TEST(SatelliteToParticleConverterTest, SatelliteToParticle) {
  // Create a particle.
  const Satellite sat{"name", SatType::DEBRIS, {10., 10., 10.}};
  const auto particle =
      LADDS::SatelliteToParticleConverter::convertSatelliteToParticle(sat, 2.2, std::numeric_limits<size_t>::max());

  const auto &expectedPosition = autopas::utils::ArrayMath::mulScalar(sat.getPosition(), 1. / 1000.0);
  const auto &expectedVelocity = autopas::utils::ArrayMath::mulScalar(sat.getVelocity(), 1. / 1000.0);

  EXPECT_THAT(particle.getPosition(), ::testing::ElementsAreArray(expectedPosition));
  EXPECT_THAT(particle.getVelocity(), ::testing::ElementsAreArray(expectedVelocity));
  EXPECT_EQ(particle.getMass(), sat.getMass());

  const auto satConverted = LADDS::SatelliteToParticleConverter::convertParticleToSatellite(particle);
  EXPECT_THAT(satConverted.getPosition(), ::testing::ElementsAreArray(sat.getPosition()));
  EXPECT_THAT(satConverted.getVelocity(), ::testing::ElementsAreArray(sat.getVelocity()));
  EXPECT_EQ(satConverted.getMass(), sat.getMass());
}

TEST(SatelliteToParticleConverterTest, ParticleToSatellite) {
  // Create a particle.
  const LADDS::Particle particle{{10., 11., 12},
                                 {1., 2., 3.},
                                 42,
                                 "dummy",
                                 LADDS::Particle::ActivityState::passive,
                                 1.,
                                 1.,
                                 LADDS::Particle::calculateBcInv(0., 1., 1., 2.2),
                                 std::numeric_limits<size_t>::max()};
  const auto sat = LADDS::SatelliteToParticleConverter::convertParticleToSatellite(particle);

  const auto &expectedPosition = autopas::utils::ArrayMath::mulScalar(particle.getPosition(), 1000.0);
  const auto &expectedVelocity = autopas::utils::ArrayMath::mulScalar(particle.getVelocity(), 1000.0);

  EXPECT_THAT(sat.getPosition(), ::testing::ElementsAreArray(expectedPosition));
  EXPECT_THAT(sat.getVelocity(), ::testing::ElementsAreArray(expectedVelocity));
  EXPECT_EQ(sat.getMass(), particle.getMass());

  const auto particleConverted =
      LADDS::SatelliteToParticleConverter::convertSatelliteToParticle(sat, 2.2, std::numeric_limits<size_t>::max());
  EXPECT_THAT(particleConverted.getPosition(), ::testing::ElementsAreArray(particle.getPosition()));
  EXPECT_THAT(particleConverted.getVelocity(), ::testing::ElementsAreArray(particle.getVelocity()));
  EXPECT_EQ(particleConverted.getMass(), particle.getMass());
}
