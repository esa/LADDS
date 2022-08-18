/**
 * @file DistanceApproximation.h
 * @author P. Gomez
 * @date 12.08.22
 */

#pragma once

namespace LADDS::utils {

/**
 * Computes distance between two particles in some interval using a linear interpolation.
 * @param i First particle
 * @param j Second particle
 * @param dt Time over which to compute the minimal distance
 * @return Vector between particles and position of p2 at minimal distance
 */
inline std::tuple<std::array<double, 3>, std::array<double, 3>> distanceByLinearInterpolation(const Particle &i,
                                                                                              const Particle &j,
                                                                                              double dt)

}  // namespace LADDS