#include "SatelliteToParticleConverterTest.h"

#include <gmock/gmock-matchers.h>

#include "ladds/SatelliteToParticleConverter.h"

#include <breakupModel/model/Satellite.h>
//#include <breakupModel/model/SatelliteBuilder.h>
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

    EXPECT_THAT(particle.getPosition(), ::testing::UnorderedElementsAreArray(expectedPosition));
    EXPECT_THAT(particle.getVelocity(), ::testing::UnorderedElementsAreArray(expectedVelocity));
}

