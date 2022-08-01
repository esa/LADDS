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
   * @param config: The configuration for defining the decomposition properties
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
   * Send the given list of leaving particles to the respective target ranks and receive their leaving
   * particles which are relevant for the local rank.
   * @param leavingParticles in/out parameter of leaving particles. If everything worked the vector should be empty
   * after the function call.
   * @param autopas
   * @param decomposition
   * @return Vector of incoming particles.
   */
  std::vector<LADDS::Particle> communicateParticles(std::vector<LADDS::Particle> &leavingParticles,
                                                    autopas::AutoPas<Particle> &autopas,
                                                    const DomainDecomposition &decomposition) const override;
};
}  // namespace LADDS