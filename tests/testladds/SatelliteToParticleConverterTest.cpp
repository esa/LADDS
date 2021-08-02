#include "SatelliteToParticleConverterTest.h"

#include <gmock/gmock-matchers.h>

#include "ladds/SatelliteToParticleConverter.h"

#include <breakupModel/model/Satellite.h>
/**
 * Sets up a particle and tries to convert it to a satellites.
 * 
 */
TEST(SatelliteToParticleConverterTest, SatelliteToParticle) {
    // Create a particle.
    const Satellite sat{"name", SatType::DEBRIS, {10.,10.,10.}};
    const auto particle = SatelliteToParticleConverter::convertSatelliteToParticle(sat);

    const auto& expectedPosition = autopas::utils::ArrayMath::mulScalar(sat.getPosition(), 1./1000.0);
    const auto& expectedVelocity = autopas::utils::ArrayMath::mulScalar(sat.getVelocity(), 1./1000.0);

    EXPECT_THAT(particle.getPosition(), ::testing::ElementsAreArray(expectedPosition));
    EXPECT_THAT(particle.getVelocity(), ::testing::ElementsAreArray(expectedVelocity));

    const auto satConverted = SatelliteToParticleConverter::convertParticleToSatellite(particle);
    EXPECT_THAT(satConverted.getPosition(), ::testing::ElementsAreArray(sat.getPosition()));
    EXPECT_THAT(satConverted.getVelocity(), ::testing::ElementsAreArray(sat.getVelocity()));
}

TEST(SatelliteToParticleConverterTest, ParticleToSatellite) {
  // Create a particle.
  const Particle particle{{10.,11.,12}, {1.,2.,3.}, 42};
  const auto sat = SatelliteToParticleConverter::convertParticleToSatellite(particle);

  const auto& expectedPosition = autopas::utils::ArrayMath::mulScalar(particle.getPosition(), 1000.0);
  const auto& expectedVelocity = autopas::utils::ArrayMath::mulScalar(particle.getVelocity(), 1000.0);

  EXPECT_THAT(sat.getPosition(), ::testing::ElementsAreArray(expectedPosition));
  EXPECT_THAT(sat.getVelocity(), ::testing::ElementsAreArray(expectedVelocity));

  const auto particleConverted = SatelliteToParticleConverter::convertSatelliteToParticle(sat);
  EXPECT_THAT(particleConverted.getPosition(), ::testing::ElementsAreArray(particle.getPosition()));
  EXPECT_THAT(particleConverted.getVelocity(), ::testing::ElementsAreArray(particle.getVelocity()));
}
