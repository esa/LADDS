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
   * Constructor.
   * @param config: The configuration for definig the decomposition properties
   */
  explicit AltitudeBasedDecomposition(ConfigReader &config);

  /**
   * Destructor.
   */
  ~AltitudeBasedDecomposition() override = default;

  /**
   * Get the rank id that contains the given simulation domain coordinates
   * @param coordinates 3D
   * @return rank id
   */
  int getRank(const std::array<double, 3> &coordinates) const override;

  /**
   * Get the particles which are leaving the local domain.
   * @param  autopas autopas container
   * @retval vector of particles
   */
  std::vector<Particle> getAndRemoveLeavingParticles(AutoPas_t &autopas) const override;

  /**
   * Get information about the grid structure
   * @return tuple{dimensions, periods, grid coordinates of current rank}
   */
  std::tuple<std::array<int, 1>, std::array<int, 1>, std::array<int, 1>> getGridInfo() const;

  /**
   * Altitude bucket boundaries for all ranks.
   */
  std::vector<double> altitudeIntervals;

  /**
   * Returns the altitude of the rank with the given index.
   * @param rank: The rank index
   * @return The altitude of the rank with the given index.
   */
  double getAltitudeOfRank(const int rank) const;

  /**
   * Creates a vector of logspaced numbers similar to np.logspace.
   * @retval std::vector<double>: The vector of logspaced numbers.
   */
  [[nodiscard]] std::vector<double> logspace(const double a, const double b, const int k);
};
}  // namespace LADDS