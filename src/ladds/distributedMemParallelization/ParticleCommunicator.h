/**
 * @file ParticleCommunicator.h
 * @author F. Gratl
 * @date 01.06.22
 */

#pragma once

#include <autopas/utils/WrapMPI.h>

#include <vector>

#include "ladds/particle/Particle.h"

namespace LADDS {

class ParticleCommunicator {
 private:
  /**
   * Vector of pairs of requests and their corresponding buffers.
   * They need to be kept in memory until their respective send operation is complete.
   */
  std::vector<std::pair<autopas::AutoPas_MPI_Request, std::vector<char>>> sendBuffers{};

 public:
  /**
   * Send a vector of particles
   * @param particles
   * @param receiver
   * @param communicator
   * @return
   */
  autopas::AutoPas_MPI_Request sendParticles(const std::vector<Particle>::iterator &particlesBegin,
                                             const std::vector<Particle>::iterator &particlesEnd,
                                             int receiver,
                                             const autopas::AutoPas_MPI_Comm &communicator);

  /**
   * Receive a vector of particles
   * @param sender
   * @param communicator
   * @return
   */
  std::vector<Particle> receiveParticles(int sender, const autopas::AutoPas_MPI_Comm &communicator);

  /**
   * Wait for all send operations to complete and then reset the buffers.
   */
  void waitAndFlushBuffers();
};
}  // namespace LADDS