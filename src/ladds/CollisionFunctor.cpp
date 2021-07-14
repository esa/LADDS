/**
 * @file CollisionFunctor.cpp
 * @author F. Gratl
 * @date 28.06.21
 */

#include "CollisionFunctor.h"

CollisionFunctor::CollisionFunctor(double cutoff) : Functor(cutoff), _cutoffSquare(cutoff * cutoff) {}

const std::unordered_map<Debris *, Debris *> &CollisionFunctor::getCollisions() const { return _collisions; }

void CollisionFunctor::AoSFunctor(Debris &i, Debris &j, bool newton3) {
  // skip interaction with deleted particles
  if (i.isDummy() or j.isDummy()) {
    return;
  }

  // calculate distance between particles
  auto dr = autopas::utils::ArrayMath::sub(i.getR(), j.getR());
  auto distanceSquare = autopas::utils::ArrayMath::dot(dr, dr);

  // if distance is wider than cutoff everything is fine and we can return
  if (distanceSquare > _cutoffSquare) {
    return;
  }

  // store pointers to colliding pair
  if (i.getID() < j.getID()) {
    _collisions[i.get<Debris::AttributeNames::ptr>()] = j.get<Debris::AttributeNames::ptr>();
  } else {
    _collisions[j.get<Debris::AttributeNames::ptr>()] = i.get<Debris::AttributeNames::ptr>();
  }
}

void CollisionFunctor::SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) {
  SoAFunctorPair(soa, soa, newton3);
}

void CollisionFunctor::SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                                      bool newton3) {
  if (soa1.getNumParticles() == 0 or soa2.getNumParticles() == 0) return;

  // get pointers to the SoA
  const auto *const __restrict ownedStatePtr1 = soa1.template begin<Debris::AttributeNames::ownershipState>();

  // outer loop over SoA1
  for (size_t i = 0; i < soa1.getNumParticles(); ++i) {
    if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
      // If the i-th particle is a dummy, skip this loop iteration.
      continue;
    }

    // inner loop over SoA2
    // TODO this one should be vectorized
#pragma omp declare reduction(vecMerge : std::vector<std::pair<Debris *, Debris *>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp declare reduction(mapMerge : std::unordered_map<Debris *, Debris *> : omp_out.insert(omp_in.begin(), omp_in.end()))
#pragma omp simd reduction(mapMerge : _collisions)
    for (size_t j = 0; j < soa2.getNumParticles(); ++j) {
      SoAKernel(i, j, soa1, soa2, newton3);
    }
  }
}

void CollisionFunctor::SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                                        bool newton3) {
  const auto *const __restrict ownedStatePtr = soa.template begin<Debris::AttributeNames::ownershipState>();

  if (ownedStatePtr[indexFirst] == autopas::OwnershipState::dummy) {
    return;
  }

  // inner loop. Should be vectorized
  for (const auto j : neighborList) {
    SoAKernel(indexFirst, j, soa, soa, newton3);
  }
}

void CollisionFunctor::SoAKernel(size_t i, size_t j, autopas::SoAView<SoAArraysType> &soa1,
                                 autopas::SoAView<SoAArraysType> &soa2, bool newton3) {
  // get pointers to the SoAs
  const auto *const __restrict x1ptr = soa1.template begin<Debris::AttributeNames::posX>();
  const auto *const __restrict y1ptr = soa1.template begin<Debris::AttributeNames::posY>();
  const auto *const __restrict z1ptr = soa1.template begin<Debris::AttributeNames::posZ>();
  const auto *const __restrict x2ptr = soa2.template begin<Debris::AttributeNames::posX>();
  const auto *const __restrict y2ptr = soa2.template begin<Debris::AttributeNames::posY>();
  const auto *const __restrict z2ptr = soa2.template begin<Debris::AttributeNames::posZ>();

  const auto *const __restrict ownedStatePtr2 = soa2.template begin<Debris::AttributeNames::ownershipState>();

  const auto *const __restrict ptr1ptr = soa1.template begin<Debris::AttributeNames::ptr>();
  const auto *const __restrict ptr2ptr = soa2.template begin<Debris::AttributeNames::ptr>();

  const auto *const __restrict id1ptr = soa1.template begin<Debris::AttributeNames::id>();
  const auto *const __restrict id2ptr = soa2.template begin<Debris::AttributeNames::id>();

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
  if (id1ptr[i] < id2ptr[j]) {
    _collisions[ptr1ptr[i]] = ptr2ptr[j];
  } else {
    _collisions[ptr2ptr[j]] = ptr1ptr[i];
  }
}
