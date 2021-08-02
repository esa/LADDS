#pragma once

#include <autopas/utils/ArrayMath.h>
#include <breakupModel/input/TLESatcatDataReader.h>
#include <breakupModel/model/Satellite.h>


#include "Particle.h"

namespace SatelliteToParticleConverter {

    /**
     * This function is used to convert the satellite data into particles.
     * @param  satellite: Satellite to convert.
     * @retval Converted particles.
     */
    [[nodiscard]] inline Particle convertSatelliteToParticle(const Satellite& satellite) {
        //Convert all entries from meters to kilometers
        const auto& position = autopas::utils::ArrayMath::mulScalar(satellite.getPosition(), 1./1000.0);
        const auto& velocity = autopas::utils::ArrayMath::mulScalar(satellite.getVelocity(), 1./1000.0);

        return Particle(position,velocity, satellite.getId());
    }

    /**
     * This function is used to convert the particle data into satellites.
     * @param  particles: Particles to convert.
     * @retval Converted satellites.
     */
    [[nodiscard]] inline Satellite convertParticleToSatellite(const Particle& particle) {

      SatelliteBuilder satelliteBuilder;

      satelliteBuilder.setID(particle.getID())
                      .setPosition(autopas::utils::ArrayMath::mulScalar(particle.getPosition(), 1000.0))
                      .setVelocity(autopas::utils::ArrayMath::mulScalar(particle.getVelocity(), 1000.0))
                      .setMass(1);

      return satelliteBuilder.getResult();
    }
}