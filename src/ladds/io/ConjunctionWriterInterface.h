/**
 * @file ConjunctionWriterInterface.h
 * @author F. Gratl
 * @date 20.12.21
 */

#pragma once

#include "ladds/particle/Particle.h"
#include "unordered_map"

class ConjuctionWriterInterface {
 public:
  ConjuctionWriterInterface() = default;

  virtual ~ConjuctionWriterInterface() = default;
  /**
   * Log the conjuction data.
   * @param iteration
   * @param collisions
   */
  virtual void writeConjunctions(size_t iteration,
                                 const std::unordered_map<Particle *, std::tuple<Particle *, double>> &collisions) = 0;
};