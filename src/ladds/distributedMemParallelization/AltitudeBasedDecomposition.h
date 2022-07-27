/**
 * @file AltitudeBasedDecomposition.h
 * @author P. Gomez
 * @date 14.07.22
 */

#pragma once

#include "DomainDecomposition.h"
#include "ladds/io/ConfigReader.h"
#include "ladds/particle/Particle.h"

namespace LADDS {

class AltitudeBasedDecomposition : public DomainDecomposition {
 public:
  /**
   * Constructor. Creates an altitude-based decomposition of particles.
   * Initially uses log-spaced intervals. A call to rebalanceDecomposition is recommended after
   * to balance the decomposition such that the number of particles per rank is approximately
   * equal.
   * @param config: The configuration for defining the decomposition properties.
   */
  explicit AltitudeBasedDecomposition(ConfigReader &config);

  /**
   * Destructor.
   */
  ~AltitudeBasedDecomposition() override = default;

  /**
   * Get the rank id that contains the given simulation domain coordinates.
   * @param coordinates 3D
   * @return rank id
   */
  int getRank(const std::array<double, 3> &coordinates) const override;

  /**
   * Get the particles that are leaving the local domain and remove them from the container.
   * These particles are not necessarily leaving the bounding box of the AutoPas container
   * but the altitude region it is responsible for.
   * @param  autopas autopas container
   * @return Vector of particles.
   */
  std::vector<Particle> getAndRemoveLeavingParticles(AutoPas_t &autopas) const override;

  /**
   * Get information about the this rank
   * @return index of current rank
   */
  int getRank() const;

  /**
   * Altitude bucket boundaries for all ranks.
   */
  std::vector<double> altitudeIntervals{};

  /**
   * Returns the altitude of the rank with the given index.
   * @param rank: The rank index
   * @return The altitude of the rank with the given index.
   */
  double getAltitudeOfRank(int rank) const;

  /**
   * Creates a vector of logspaced numbers similar to np.logspace.
   * @return Vector of logspaced numbers.
   */
  [[nodiscard]] std::vector<double> logspace(const double a, const double b, const int k);

  /**
   * Balances the domain decomposition based on particle locations,
   * aiming for a similar number of particles per rank.
   * @param particles vector of all particles in all ranks
   * @param autopas
   */
  void rebalanceDecomposition(const std::vector<LADDS::Particle> &particles, AutoPas_t &autopas) override;
};
}  // namespace LADDS