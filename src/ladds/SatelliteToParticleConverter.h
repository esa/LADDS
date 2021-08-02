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
        throw std::runtime_error("Not implemented. Ask Future Fabio about this.");
    }
}