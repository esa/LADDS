/**
 * @file RegularGridDecomposition.h
 * @author F. Gratl
 * @date 20.05.22
 */

#pragma once

#include "DomainDecomposition.h"
#include "ladds/io/ConfigReader.h"

namespace LADDS {
class RegularGridDecomposition : public DomainDecomposition {
 public:
  /**
   * Constructor.
   * @param config: The configuration for definig the decomposition properties
   */
  explicit RegularGridDecomposition(ConfigReader &config);

  /**
   * Destructor.
   */
  ~RegularGridDecomposition() override = default;

  int getRank(const std::array<double, 3> &coordinates) const override;

  /**
   * Get information about the grid structure
   * @return tuple{dimensions, periods, grid coordinates of current rank}
   */
  std::tuple<std::array<int, 3>, std::array<int, 3>, std::array<int, 3>> getGridInfo() const;

  /**
   * Get the particles which are leaving the local domain.
   * @param  autopas autopas container
   * @retval vector of particles
   */
  std::vector<Particle> getLeavingParticles(const AutoPas_t &autopas) const override;
};
}  // namespace LADDS