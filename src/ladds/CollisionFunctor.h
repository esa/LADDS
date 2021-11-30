/**
 * @file CollisionFunctor.h
 * @author F. Gratl
 * @date 28.06.21
 */

#pragma once

#include <autopas/pairwiseFunctors/Functor.h>
#include <autopas/utils/SoA.h>
#include <autopas/utils/SoAView.h>

#include "ladds/particle/Particle.h"

/**
 * Class describing the pairwise particle computation for determining whether a collision occured.
 * This functor is passed to autopas::AutoPas::iteratePairwise() as the primary pairwise interaction.
 */
class CollisionFunctor final : public autopas::Functor<Particle, CollisionFunctor> {
 public:
  explicit CollisionFunctor(double cutoff);

  [[nodiscard]] bool isRelevantForTuning() final {
    return true;
  }

  [[nodiscard]] bool allowsNewton3() final {
    return true;
  }

  [[nodiscard]] bool allowsNonNewton3() final {
    return true;
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

  [[nodiscard]] const std::unordered_map<Particle *, Particle *> &getCollisions() const;

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
    std::unordered_map<Particle *, Particle *> collisions{};
  } __attribute__((aligned(64)));

  // make sure that the size of ThreadData is correct
  static_assert(sizeof(ThreadData) % 64 == 0, "ThreadData has wrong size");

  // Buffer per thread
  std::vector<ThreadData> _threadData{};
  // key = particle with the smaller id
  std::unordered_map<Particle *, Particle *> _collisions{};
  const double _cutoffSquare;
};
