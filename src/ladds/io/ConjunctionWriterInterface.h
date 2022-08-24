/**
 * @file ConjunctionWriterInterface.h
 * @author F. Gratl
 * @date 20.12.21
 */

#pragma once

#include "ladds/CollisionFunctor.h"
#include "ladds/particle/Particle.h"

namespace LADDS {

class ConjuctionWriterInterface {
 public:
  ConjuctionWriterInterface() = default;

  virtual ~ConjuctionWriterInterface() = default;
  /**
   * Log the conjuction data.
   * @param iteration
   * @param collisions
   */
  virtual void writeConjunctions(size_t iteration, const CollisionFunctor::CollisionCollectionT &collisions) = 0;

  /**
   * Log the conjuction data.
   * @param iteration
   * @param collisions
   */
  virtual void writeEvasions(size_t iteration, const CollisionFunctor::CollisionCollectionT &evasions) = 0;
};

}  // namespace LADDS