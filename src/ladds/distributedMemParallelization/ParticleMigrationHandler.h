/**
 * @file ParticleMigrationHandler.h
 * @author P. Gomez
 * @date 18.07.22
 */
#pragma once

#include <autopas/AutoPasDecl.h>

#include <vector>

#include "ladds/CollisionFunctor.h"
#include "ladds/TypeDefinitions.h"
#include "ladds/distributedMemParallelization/DomainDecomposition.h"
#include "ladds/particle/Particle.h"

/**
 * Provides tools for sending particles between ranks and
 * detecting conjunctions at the boundaries of the domain.
 */
namespace LADDS::ParticleMigrationHandler {

/**
 * Interact all incoming particles with all particles which potentially crossed its path since the last container
 * update.
 *
 * These particles are found in a box around the immigrant's position -deltaT time ago.
 * The box has a side length of 2x the maximum coverable distance by any particle.
 *
 * @note The first pointers in the returned tuple collection point to particles in the immigrant vector!
 *
 * @param autopas
 * @param incomingParticles
 * @param deltaT Time since the last container update.
 * @param maxV Maximal velocity a particle is assumed to have. Has to be positive.
 * @param collisionDistanceFactor See CollisionFunctor::_collisionDistanceFactor
 * @param minDetectionRadius See CollisionFunctor::_minDetectionRadius
 * @param CDMcutoffInKM See CollisionFunctor::_CDMcutoffInKM
 * @return Collection of collision partners
 */
std::tuple<LADDS::CollisionFunctor::CollisionCollectionT, LADDS::CollisionFunctor::CollisionCollectionT>
collisionDetectionImmigrants(AutoPas_t &autopas,
                             std::vector<LADDS::Particle> &incomingParticles,
                             double deltaT,
                             double maxV,
                             double collisionDistanceFactor,
                             double minDetectionRadius,
                             double CDMcutoffInKM);

}  // namespace LADDS::ParticleMigrationHandler