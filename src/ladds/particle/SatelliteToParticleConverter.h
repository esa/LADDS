#pragma once

#include <autopas/utils/ArrayMath.h>
#include <breakupModel/input/TLESatcatDataReader.h>
#include <breakupModel/model/Satellite.h>

#include "Particle.h"

namespace SatelliteToParticleConverter {

/**
 * This function is used to convert the satellite data into particles.
 * @note Resulting particles always have Particle::ActivityState::passive as we consider results of a collision to be
 * broken. If a non-catastrophic collision should be simulated the Particle::ActivityState has to be changed on the
 * returned object.
 * @param  satellite: Satellite to convert.
 * @retval Converted particles.
 */
[[nodiscard]] inline Particle convertSatelliteToParticle(const Satellite &satellite) {
  // Convert all entries from meters to kilometers
  const auto &position = autopas::utils::ArrayMath::mulScalar(satellite.getPosition(), 1. / 1000.0);
  const auto &velocity = autopas::utils::ArrayMath::mulScalar(satellite.getVelocity(), 1. / 1000.0);
  const auto radius = std::sqrt(satellite.getArea() * M_1_PI);
  return Particle{
      position, velocity, satellite.getId(), Particle::ActivityState::passive, satellite.getMass(), radius, 2.2};
}

/**
 * This function is used to convert the particle data into satellites.
 * @note This conversion does not retain the information about the satellites radius and activity state!
 * @param  particles: Particles to convert.
 * @retval Converted satellites.
 */
[[nodiscard]] inline Satellite convertParticleToSatellite(const Particle &particle) {
  SatelliteBuilder satelliteBuilder;

  satelliteBuilder
      .setID(particle.getID())
      // Converting from km to meters
      .setPosition(autopas::utils::ArrayMath::mulScalar(particle.getPosition(), 1000.0))
      .setVelocity(autopas::utils::ArrayMath::mulScalar(particle.getVelocity(), 1000.0))
      .setMass(particle.getMass());

  return satelliteBuilder.getResult();
}
}  // namespace SatelliteToParticleConverter
