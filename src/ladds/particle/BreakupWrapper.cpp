/**
 * @file BreakupWrapper.cpp
 * @author F. Gratl
 * @date 24.01.22
 */

#include "BreakupWrapper.h"

void BreakupWrapper::simulateBreakup(const CollisionFunctor::CollisionCollectionT &collisions) {
  for (const auto &[p1, p2, dist] : collisions) {
    // convert particles to satellites for breakup code
    std::vector<Satellite> satellites{
        SatelliteToParticleConverter::convertParticleToSatellite(*p1),
        SatelliteToParticleConverter::convertParticleToSatellite(*p2),
    };
    // remove colliding particles from autopas. If the collision is non-fatal they will be reinserted further down
    autopas.deleteParticle(*p1);
    autopas.deleteParticle(*p2);

    // Create and run the Breakup Simulation
    auto configSource = std::make_shared<RuntimeInputSource>(
        minLc, satellites, SimulationType::COLLISION, maxExistingParticleId, std::nullopt, true);
    BreakupBuilder breakupBuilder{configSource};
    auto breakup = breakupBuilder.getBreakup();
    breakup->run();

    // insert resulting debris into the simulation
    const auto &newSatellites = breakup->getResult();
    for (const auto &newSat : newSatellites) {
      autopas.addParticle(SatelliteToParticleConverter::convertSatelliteToParticle(newSat));
    }
  }
}
