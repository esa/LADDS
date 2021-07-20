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

  [[nodiscard]] static constexpr auto getNeededAttr() {
    return std::array<typename Debris::AttributeNames, 6>{
        Debris::AttributeNames::ptr,  Debris::AttributeNames::id,   Debris::AttributeNames::ownershipState,
        Debris::AttributeNames::posX, Debris::AttributeNames::posY, Debris::AttributeNames::posZ};
  }

  [[nodiscard]] static constexpr auto getNeededAttr(std::false_type) { return getNeededAttr(); }

  [[nodiscard]] static constexpr std::array<typename Debris::AttributeNames, 0> getComputedAttr() {
    return std::array<typename Debris::AttributeNames, 0>{/*Nothing*/};
  };

  void initTraversal() final;

  void endTraversal(bool newton3) final;

  [[nodiscard]] const std::unordered_map<Debris *, Debris *> &getCollisions() const;

  void AoSFunctor(Debris &i, Debris &j, bool newton3) final;

  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final;

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2, bool newton3) final;

  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList, bool newton3) final;

 private:
  void SoAKernel(size_t i, size_t j, autopas::SoAView<SoAArraysType> &soa1, autopas::SoAView<SoAArraysType> &soa2,
                 bool newton3);

  // Buffer struct that is safe against false sharing
  //! @cond Doxygen_Suppress
  struct ThreadData {
    std::unordered_map<Debris *, Debris *> collisions{};
  } __attribute__((aligned(64)));

  // make sure that the size of ThreadData is correct
  static_assert(sizeof(ThreadData) % 64 == 0, "ThreadData has wrong size");
  //! @endcond

  // Buffer per thread
  std::vector<ThreadData> _threadData;
  // key = particle with the smaller id
  std::unordered_map<Debris *, Debris *> _collisions{};
  const double _cutoffSquare;
};
