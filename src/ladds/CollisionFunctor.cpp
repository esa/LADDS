/**
 * @file CollisionFunctor.cpp
 * @author F. Gratl
 * @date 28.06.21
 */

#include "CollisionFunctor.h"

#include <autopas/utils/ArrayMath.h>
#include <autopas/utils/WrapOpenMP.h>

#include <algorithm>

CollisionFunctor::CollisionFunctor(double cutoff, double dt, double minorCutoff)
    : Functor(cutoff), _cutoffSquare(cutoff * cutoff), _dt(dt), _minorCutoffSquare(minorCutoff * minorCutoff) {
  _threadData.resize(autopas::autopas_get_max_threads());
}

const CollisionFunctor::CollisionCollectionT &CollisionFunctor::getCollisions() const {
  return _collisions;
}

void CollisionFunctor::AoSFunctor(Particle &i, Particle &j, bool newton3) {
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::dot;
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::sub;

  // skip interaction with deleted particles
  if (i.isDummy() or j.isDummy()) {
    return;
  }

  // calculate distance between particles
  const auto dr = sub(i.getR(), j.getR());
  const auto distanceSquare = dot(dr, dr);

  // if distance is wider than cutoff everything is fine and we can return
  if (distanceSquare > _cutoffSquare) {
    return;
  }

  // if distance is smaller we compute the minimum distance between the linear interpolations
  // of the particle trajectory over the last timestep
  // according to wolfram alpha, should look like this:
  // https://www.wolframalpha.com/input/?i=solve+for+t+d%2Fdt%5Bsqrt%28%28-z_1+t+%2B+t+v_1+%2B+x_1+-+y_1%29%5E2+%2B+%28-z_2+t+%2B+t+v_2+%2B+x_2+-+y_2%29%5E2+%2B+%28t+v_3+-+t+z_3+%2B+x_3+-+y_3%29%5E2%29%5D

  const auto &vi = i.getVelocity();
  const auto &vj = j.getVelocity();

  // Get old time step position
  const auto old_r_i = sub(i.getR(), mulScalar(vi, _dt));
  const auto old_r_j = sub(j.getR(), mulScalar(vj, _dt));

  // Compute nominator dot products
  const auto vi_ri = dot(vi, old_r_i);
  const auto vi_rj = dot(vi, old_r_j);
  const auto vj_ri = dot(vj, old_r_i);
  const auto vj_rj = dot(vj, old_r_j);

  const auto nominator = vi_rj + vj_ri - vi_ri - vj_rj;

  // Compute denominator dot products
  const auto two_vi_vj = 2.0 * dot(vi, vj);
  const auto vi_square = dot(vi, vi);
  const auto vj_square = dot(vj, vj);

  const auto denominator = vi_square + vj_square - two_vi_vj;

  // Compute t at minimal distance
  auto t = nominator / denominator;

  // If in the past, minimum is at t = 0
  // Else If in future timesteps, minimum for this is at t = _dt
  t = std::clamp(t, 0., _dt);

  // Compute actual distance by propagating along the line to t
  const auto p1 = add(old_r_i, mulScalar(vi, t));
  const auto p2 = add(old_r_j, mulScalar(vj, t));

  const auto dr_lines = sub(p1, p2);
  const auto distanceSquare_lines = dot(dr_lines, dr_lines);

  if (distanceSquare_lines > _minorCutoffSquare) return;

  // store pointers to colliding pair
  if (i.getID() < j.getID()) {
    _threadData[autopas::autopas_get_thread_num()].collisions.emplace_back(
        i.get<Particle::AttributeNames::ptr>(), j.get<Particle::AttributeNames::ptr>(), distanceSquare_lines);
  } else {
    _threadData[autopas::autopas_get_thread_num()].collisions.emplace_back(
        j.get<Particle::AttributeNames::ptr>(), i.get<Particle::AttributeNames::ptr>(), distanceSquare_lines);
  }
}

void CollisionFunctor::SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) {
  SoAFunctorPair(soa, soa, newton3);
}

void CollisionFunctor::SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1,
                                      autopas::SoAView<SoAArraysType> soa2,
                                      bool newton3) {
  if (soa1.getNumParticles() == 0 or soa2.getNumParticles() == 0) return;

  // get pointers to the SoA
  const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle::AttributeNames::ownershipState>();

  // outer loop over SoA1
  for (size_t i = 0; i < soa1.getNumParticles(); ++i) {
    if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
      // If the i-th particle is a dummy, skip this loop iteration.
      continue;
    }

    // inner loop over SoA2
    // custom reduction for collision collections (vector)
#pragma omp declare reduction(vecMerge:CollisionCollectionT \
                              : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
    // alias because OpenMP needs it
    auto &thisCollisions = _threadData[autopas::autopas_get_thread_num()].collisions;
#pragma omp simd reduction(vecMerge : thisCollisions)
    for (size_t j = 0; j < soa2.getNumParticles(); ++j) {
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
      "SoA kernel not up to date with AoS Kernel as it lacks the new linear interpolation distance.");

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
  if (id1ptr[i] < id2ptr[j]) {
    _threadData[autopas::autopas_get_thread_num()].collisions.emplace_back(ptr1ptr[i], ptr2ptr[j], dr2);
  } else {
    _threadData[autopas::autopas_get_thread_num()].collisions.emplace_back(ptr2ptr[j], ptr1ptr[i], dr2);
  }
}
void CollisionFunctor::initTraversal() {
  _collisions.clear();
}

void CollisionFunctor::endTraversal(bool newton3) {
  for (auto &data : _threadData) {
    _collisions.insert(_collisions.end(), data.collisions.begin(), data.collisions.end());
    data.collisions.clear();
  }
}
