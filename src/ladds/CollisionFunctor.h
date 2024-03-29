/**
 * @file CollisionFunctor.h
 * @author F. Gratl
 * @date 28.06.21
 */

#pragma once

#include <autopas/pairwiseFunctors/Functor.h>
#include <autopas/utils/SoA.h>
#include <autopas/utils/SoAView.h>

#include <array>

#include "ladds/particle/Particle.h"

namespace LADDS {
/**
 * Class describing the pairwise particle computation for determining whether a collision occured.
 * This functor is passed to autopas::AutoPas::iteratePairwise() as the primary pairwise interaction.
 */
class CollisionFunctor final : public autopas::Functor<Particle, CollisionFunctor> {
 public:
  /**
   * Constructor
   * @param cutoff Distance for two particles to be considered for sub time step investigation (interpolation).
   * @param dt time over which to consider collisions (length of the interpolation).
   * @param collisionDistanceFactor See CollisionFunctor::_collisionDistanceFactor.
   * @param minDetectionRadius All particles with a larger radius are assumed to be detectable by radar.
   *        Thus collisions with particles that are Particle::ActivityState::evasive will not be considered.
   * @param evasionTrackingCutoffInKM Absolute cutoff we consider for CDMs , that is tracking of evaded conjunctions.
   */
  CollisionFunctor(double cutoff,
                   double dt,
                   double collisionDistanceFactor,
                   double minDetectionRadius,
                   double evasionTrackingCutoffInKM);

  [[nodiscard]] bool isRelevantForTuning() final {
    return true;
  }

  [[nodiscard]] bool allowsNewton3() final {
    return true;
  }

  [[nodiscard]] bool allowsNonNewton3() final {
    // as we are only interested in any interaction [i,j] it makes no sense to also look at [j,i]
    return false;
  }

  [[nodiscard]] constexpr static auto getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 6>{Particle::AttributeNames::ptr,
                                                            Particle::AttributeNames::id,
                                                            Particle::AttributeNames::ownershipState,
                                                            Particle::AttributeNames::posX,
                                                            Particle::AttributeNames::posY,
                                                            Particle::AttributeNames::posZ};
  }

  [[nodiscard]] constexpr static auto getNeededAttr(std::false_type) {
    return getNeededAttr();
  }

  [[nodiscard]] constexpr static std::array<typename Particle::AttributeNames, 0> getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 0>{/*Nothing*/};
  };

  void initTraversal() final;

  void endTraversal(bool newton3) final;

  using CollisionT = std::tuple<Particle *, Particle *, double, std::array<double, 3>>;
  using CollisionCollectionT = std::vector<CollisionT>;

  [[nodiscard]] const CollisionCollectionT &getCollisions() const;

  [[nodiscard]] const CollisionCollectionT &getEvadedCollisions() const;

  void AoSFunctor(Particle &i, Particle &j, bool newton3) final;

  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final;

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2, bool newton3) final;

  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa,
                        const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final;

 private:
  void SoAKernel(
      size_t i, size_t j, autopas::SoAView<SoAArraysType> &soa1, autopas::SoAView<SoAArraysType> &soa2, bool newton3);

  // Buffer struct that is safe against false sharing
  struct ThreadData {
    CollisionCollectionT collisions{};
    CollisionCollectionT evadedCollisions{};
  } __attribute__((aligned(64)));

  // make sure that the size of ThreadData is correct
  static_assert(sizeof(ThreadData) % 64 == 0, "ThreadData has wrong size");

  // Buffer per thread
  std::vector<ThreadData> _threadData{};
  // key = particle with the smaller id
  CollisionCollectionT _collisions{};
  CollisionCollectionT _evadedCollisions{};
  const double _cutoffSquare;
  const double _dt;

  /**
   * Factor multiplied with the sum of radii to over approximate collision distances.
   * Also converts from meter to km.
   */
  const double _collisionDistanceFactor;

  /**
   * Minimal object size in meter that is assumed to be detectable.
   * Objects larger than this can be evaded by evasive particles.
   */
  const double _minDetectionRadius;

  /**
   * Squared of cutoff we consider for conjunction detection messages (CDMs) , that determines
   * whether we are tracking them as evaded conjunctions.
   */
  const double _squaredEvasionTrackingCutoffInKM;
};

}  // namespace LADDS