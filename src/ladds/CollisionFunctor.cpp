/**
 * @file CollisionFunctor.cpp
 * @author F. Gratl
 * @date 28.06.21
 */

#include "CollisionFunctor.h"

#include <autopas/utils/ArrayMath.h>
#include <autopas/utils/WrapOpenMP.h>

#include <algorithm>
#include <limits>

#include "ladds/utils/DistanceApproximation.h"

namespace LADDS {
CollisionFunctor::CollisionFunctor(
    double cutoff, double dt, double collisionDistanceFactor, double minDetectionRadius, double CDMcutoffInKM)
    : Functor(cutoff),
      _cutoffSquare(cutoff * cutoff),
      _dt(dt),
      _collisionDistanceFactor(collisionDistanceFactor / 1000.),  // also imply conversion from m to km
      _minDetectionRadius(minDetectionRadius),
      _squaredCDMcutoffInKM(CDMcutoffInKM * CDMcutoffInKM) {
  _threadData.resize(autopas::autopas_get_max_threads());
}

const CollisionFunctor::CollisionCollectionT &CollisionFunctor::getCollisions() const {
  return _collisions;
}

const CollisionFunctor::CollisionCollectionT &CollisionFunctor::getEvadedCollisions() const {
  return _evadedCollisions;
}

void CollisionFunctor::AoSFunctor(Particle &i, Particle &j, bool newton3) {
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::dot;
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::sub;

  // skip if interaction involves:
  //  - any deleted particles
  //  - two actively evasive satellites
  //  - one evasive and one of detectable size
  if (i.isDummy() or j.isDummy()) {
    return;
  }

  const auto &iActivity = i.getActivityState();
  const auto &jActivity = j.getActivityState();
  const auto &iRadius = i.getRadius();
  const auto &jRadius = j.getRadius();
  bool wasEvaded = false;
  if ((iActivity != Particle::ActivityState::passive and jActivity != Particle::ActivityState::passive) or
      (iActivity != Particle::ActivityState::passive and jRadius >= _minDetectionRadius) or
      (jActivity != Particle::ActivityState::passive and iRadius >= _minDetectionRadius)) {
    wasEvaded = true;
  }

  // skip if both particles have same parent, i.e. originate from same breakup
  if ((i.getParentIdentifier() != std::numeric_limits<size_t>::max()) and
      i.getParentIdentifier() == j.getParentIdentifier()) {
    return;
  }
  // calculate distance between particles
  const auto dr = sub(i.getR(), j.getR());
  const auto distanceSquare = dot(dr, dr);

  // if distance is wider than cutoff everything is fine and we can return
  if (distanceSquare > _cutoffSquare) {
    return;
  }

  // if distance is smaller we compute the shortest vector between the linear interpolations
  // of the particle trajectory over the last timestep and position of particle 2 at that time
  const auto [dr_lines, p2] = distanceByLinearInterpolation(i, j, _dt);

  const auto distanceSquare_lines = dot(dr_lines, dr_lines);
  // _collisionDistanceFactor also converts this to km
  const auto scaledObjectSeparation = (iRadius + jRadius) * _collisionDistanceFactor;

  // For evaded collisions we collect them at this point
  // and store pointers to colliding pair
  if (wasEvaded) {
    if (distanceSquare_lines < _squaredCDMcutoffInKM) {
      return;
    }
    // compute potential collision point as middle between the two particles
    const auto collision_point = add(p2, mulScalar(dr_lines, 0.5));

    if (i.getID() < j.getID()) {
      _threadData[autopas::autopas_get_thread_num()].evadedCollisions.emplace_back(
          i.get<Particle::AttributeNames::ptr>(),
          j.get<Particle::AttributeNames::ptr>(),
          distanceSquare_lines,
          collision_point);
    } else {
      _threadData[autopas::autopas_get_thread_num()].evadedCollisions.emplace_back(
          j.get<Particle::AttributeNames::ptr>(),
          i.get<Particle::AttributeNames::ptr>(),
          distanceSquare_lines,
          collision_point);
    }
  } else {
    // if they were not evaded, check if actual collisions takes place, i.e.
    // particles are far enough away (i.e. further than sum of their scaled sizes)
    if (distanceSquare_lines > (scaledObjectSeparation * scaledObjectSeparation)) {
      return;
    }

    // if not compute collision point as middle between the two particles
    const auto collision_point = add(p2, mulScalar(dr_lines, 0.5));

    if (i.getID() < j.getID()) {
      _threadData[autopas::autopas_get_thread_num()].collisions.emplace_back(i.get<Particle::AttributeNames::ptr>(),
                                                                             j.get<Particle::AttributeNames::ptr>(),
                                                                             distanceSquare_lines,
                                                                             collision_point);
    } else {
      _threadData[autopas::autopas_get_thread_num()].collisions.emplace_back(j.get<Particle::AttributeNames::ptr>(),
                                                                             i.get<Particle::AttributeNames::ptr>(),
                                                                             distanceSquare_lines,
                                                                             collision_point);
    }
  }
}

void CollisionFunctor::SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) {
  SoAFunctorPair(soa, soa, newton3);
}

void CollisionFunctor::SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1,
                                      autopas::SoAView<SoAArraysType> soa2,
                                      bool newton3) {
  if (soa1.getNumberOfParticles() == 0 or soa2.getNumberOfParticles() == 0) return;

  // get pointers to the SoA
  const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle::AttributeNames::ownershipState>();

  // outer loop over SoA1
  for (size_t i = 0; i < soa1.getNumberOfParticles(); ++i) {
    if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
      // If the i-th particle is a dummy, skip this loop iteration.
      continue;
    }

    // inner loop over SoA2
    // custom reduction for collision collections (vector)
#pragma omp declare reduction(vecMerge:CollisionCollectionT \
                              : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
    // alias because OpenMP needs it
    // TODO add evadedCollisions when refactoring this
    auto &thisCollisions = _threadData[autopas::autopas_get_thread_num()].collisions;
#pragma omp simd reduction(vecMerge : thisCollisions)
    for (size_t j = 0; j < soa2.getNumberOfParticles(); ++j) {
      SoAKernel(i, j, soa1, soa2, newton3);
    }
  }
}

void CollisionFunctor::SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa,
                                        const size_t indexFirst,
                                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                                        bool newton3) {
  const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

  if (ownedStatePtr[indexFirst] == autopas::OwnershipState::dummy) {
    return;
  }

  // inner loop. Should be vectorized
  for (const auto j : neighborList) {
    SoAKernel(indexFirst, j, soa, soa, newton3);
  }
}

void CollisionFunctor::SoAKernel(
    size_t i, size_t j, autopas::SoAView<SoAArraysType> &soa1, autopas::SoAView<SoAArraysType> &soa2, bool newton3) {
  // TODO: as soon as this exception is removed / the SoAFunctor properly implemented
  // remove the GTEST_SKIP in CollisionFunctorIntegrationTest::testAutoPasAlgorithm!
  throw std::runtime_error(
      "SoA kernel not up to date with AoS Kernel as it lacks the new linear interpolation distance and correct "
      "collision position.");

  // get pointers to the SoAs
  const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
  const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
  const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
  const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
  const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
  const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

  const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle::AttributeNames::ownershipState>();

  const auto *const __restrict ptr1ptr = soa1.template begin<Particle::AttributeNames::ptr>();
  const auto *const __restrict ptr2ptr = soa2.template begin<Particle::AttributeNames::ptr>();

  const auto *const __restrict id1ptr = soa1.template begin<Particle::AttributeNames::id>();
  const auto *const __restrict id2ptr = soa2.template begin<Particle::AttributeNames::id>();

  if (ownedStatePtr2[j] == autopas::OwnershipState::dummy or id1ptr[i] == id2ptr[j]) {
    return;
  }

  const auto drx = x1ptr[i] - x2ptr[j];
  const auto dry = y1ptr[i] - y2ptr[j];
  const auto drz = z1ptr[i] - z2ptr[j];

  const auto drx2 = drx * drx;
  const auto dry2 = dry * dry;
  const auto drz2 = drz * drz;

  const auto dr2 = drx2 + dry2 + drz2;

  if (dr2 > _cutoffSquare) {
    return;
  }

  // store pointers to colliding pair
  // TODO add evaded collisions when refactoring this
  if (id1ptr[i] < id2ptr[j]) {
    _threadData[autopas::autopas_get_thread_num()].collisions.emplace_back(
        ptr1ptr[i], ptr2ptr[j], dr2, ptr1ptr[i]->getPosition());
  } else {
    _threadData[autopas::autopas_get_thread_num()].collisions.emplace_back(
        ptr2ptr[j], ptr1ptr[i], dr2, ptr1ptr[i]->getPosition());
  }
}
void CollisionFunctor::initTraversal() {
  _collisions.clear();
  _evadedCollisions.clear();
}

void CollisionFunctor::endTraversal(bool newton3) {
  for (auto &data : _threadData) {
    _collisions.insert(_collisions.end(), data.collisions.begin(), data.collisions.end());
    data.collisions.clear();
    _evadedCollisions.insert(_evadedCollisions.end(), data.evadedCollisions.begin(), data.evadedCollisions.end());
    data.evadedCollisions.clear();
  }
}
}  // namespace LADDS