/**
 * @file DomainDecomposition.h
 * @author F. Gratl
 * @date 20.05.22
 */

#pragma once

#include <autopas/utils/WrapMPI.h>

#include <array>
#include <vector>

#include "ladds/TypeDefinitions.h"
#include "ladds/particle/Particle.h"

namespace LADDS {

/**
 * Interface for MPI domain decompositions.
 */
class DomainDecomposition {
 public:
  /**
   * Destructor.
   */
  virtual ~DomainDecomposition() = default;

  /**
   * Returns the minimum coordinates of the global domain.
   * @return bottom left front corner of the global domain.
   */
  [[nodiscard]] virtual std::array<double, 3> getGlobalBoxMin() const;

  /**
   * Returns the maximum coordinates of the global domain.
   * @return top right back corner of the global domain.
   */
  [[nodiscard]] virtual std::array<double, 3> getGlobalBoxMax() const;

  /**
   * Returns the minimum coordinates of the local domain.
   * @return bottom left front corner of the local domain.
   */
  [[nodiscard]] virtual std::array<double, 3> getLocalBoxMin() const;

  /**
   * Returns the maximum coordinates of the local domain.
   * @return top right back corner of the local domain.
   */
  [[nodiscard]] virtual std::array<double, 3> getLocalBoxMax() const;

  /**
   * Checks if the provided coordinates are in the global domain.
   * @param coordinates 3D
   * @return true iff the coordinates are inside the global domain.
   */
  [[nodiscard]] virtual bool containsGlobal(const std::array<double, 3> &coordinates) const;

  /**
   * Checks if the provided coordinates are in the local domain.
   * @param coordinates 3D
   * @return true iff the coordinates are inside the local domain.
   */
  [[nodiscard]] virtual bool containsLocal(const std::array<double, 3> &coordinates) const;

  /**
   * Get the rank id that contains the given simulation domain coordinates
   * @param coordinates 3D
   * @return rank id
   */
  [[nodiscard]] virtual int getRank(const std::array<double, 3> &coordinates) const = 0;

  /**
   * Get the particles which are leaving the local domain.
   * @param  autopas autopas container
   * @retval vector of particles
   */
  virtual std::vector<Particle> getAndRemoveLeavingParticles(AutoPas_t &autopas) const {
    throw std::runtime_error(
        "This function is not implemeted. Leaving particles can be identified with autopas.updateContainer");
  }

  /**
   * Get the communicator used by this decomposition.
   * @return
   */
  autopas::AutoPas_MPI_Comm getCommunicator() const;

  /**
   * Balances the domain decomposition based on particle locations.
   * @param  &particles: vector of all particles in all ranks
   */
  virtual void rebalanceDecomposition(const std::vector<LADDS::Particle> &particles, AutoPas_t &autopas) {
    // This function is currently not implemented.
  }

 protected:
  std::array<double, 3> globalBoxMin{};
  std::array<double, 3> globalBoxMax{};
  std::array<double, 3> localBoxMin{};
  std::array<double, 3> localBoxMax{};
  /**
   * The MPI communicator containing all processes which own a subdomain in this decomposition.
   */
  autopas::AutoPas_MPI_Comm communicator{};
};

}  // namespace LADDS