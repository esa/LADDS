/**
 * @file BreakupWrapper.h
 * @author F. Gratl
 * @date 24.01.22
 */

#pragma once

#include <breakupModel/input/RuntimeInputSource.h>
#include <breakupModel/simulation/BreakupBuilder.h>

#include "ladds/CollisionFunctor.h"
#include "ladds/TypeDefinitions.h"
#include "ladds/io/ConfigReader.h"
#include "ladds/particle/SatelliteToParticleConverter.h"

/**
 * Wrapper class for the NASA Breakup Model.
 * See https://github.com/esa/NASA-breakup-model-cpp
 */
class BreakupWrapper {
 public:
  BreakupWrapper(ConfigReader &config, AutoPas_t &autopas)
      : minLc{config.get<double>("sim/breakup/minLc", 0.05)},
        enforceMassConservation{config.get<bool>("sim/breakup/enforceMassConservation", true)},
        coefficientOfDrag{config.get<double>("sim/coefficientOfDrag")},
        autopas{autopas} {
    // find the highest existing particle id:
    // Particles are not sorted by id and might neither be starting by 0 nor be consecutive (e.g. due to burn-ups)
    // therefore we have to go through all of them
    for (const auto &p : autopas) {
      maxExistingParticleId = std::max(maxExistingParticleId, p.getID());
    }
  };

  void simulateBreakup(const CollisionFunctor::CollisionCollectionT &collisions);

 private:
  const double minLc;
  const bool enforceMassConservation;
  double coefficientOfDrag;
  size_t maxExistingParticleId{0ul};
  AutoPas_t &autopas;
};
