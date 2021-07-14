/**
 * @file CollisionFunctor.h
 * @author F. Gratl
 * @date 28.06.21
 */

#pragma once

#include <autopas/pairwiseFunctors/Functor.h>
#include <autopas/utils/SoA.h>
#include <autopas/utils/SoAView.h>

#include "Debris.h"

/**
 * Class describing the pairwise particle computation for determining whether a collision occured.
 * This functor is passed to autopas::AutoPas::iteratePairwise() as the primary pairwise interaction.
 */
class CollisionFunctor final : public autopas::Functor<Debris, CollisionFunctor> {
 public:
  explicit CollisionFunctor(double cutoff);

  [[nodiscard]] bool isRelevantForTuning() final { return true; }

  [[nodiscard]] bool allowsNewton3() final { return true; }

  [[nodiscard]] bool allowsNonNewton3() final { return true; }

  [[nodiscard]] constexpr static auto getNeededAttr() {
    return std::array<typename Debris::AttributeNames, 6>{
        Debris::AttributeNames::ptr,  Debris::AttributeNames::id,   Debris::AttributeNames::ownershipState,
        Debris::AttributeNames::posX, Debris::AttributeNames::posY, Debris::AttributeNames::posZ};
  }

  [[nodiscard]] constexpr static auto getNeededAttr(std::false_type) { return getNeededAttr(); }

  [[nodiscard]] constexpr static std::array<typename Debris::AttributeNames, 0> getComputedAttr() {
    return std::array<typename Debris::AttributeNames, 0>{/*Nothing*/};
  };

  void initTraversal() final{};

  void endTraversal(bool newton3) final{};

  [[nodiscard]] const std::unordered_map<Debris *, Debris *> &getCollisions() const;

  void AoSFunctor(Debris &i, Debris &j, bool newton3) final;

  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final;

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2, bool newton3) final;

  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList, bool newton3) final;

 private:
  void SoAKernel(size_t i, size_t j, autopas::SoAView<SoAArraysType> &soa1, autopas::SoAView<SoAArraysType> &soa2,
                 bool newton3);

  // TODO make this thread-safe, false sharing, etc...
  // key = particle with the smaller id
  std::unordered_map<Debris *, Debris *> _collisions{};
  const double _cutoffSquare;
};