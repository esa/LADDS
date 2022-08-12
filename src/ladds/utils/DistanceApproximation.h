/**
 * @file DistanceApproximation.h
 * @author P. Gomez
 * @date 12.08.22
 */

#pragma once

#include <autopas/utils/ArrayMath.h>

namespace LADDS {

/**
 * Computes distance between two particles in some interval using a linear interpolation.
 * @param i First particle
 * @param j Second particle
 * @param dt Time over which to compute the minimal distance
 * @return Vector between particles and position of p2 at minimal distance
 */
inline std::tuple<std::array<double, 3>, std::array<double, 3>> distanceByLinearInterpolation(const Particle &i,
                                                                                              const Particle &j,
                                                                                              double dt) {
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::dot;
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::sub;

  // according to wolfram alpha, should look like this:
  // https://www.wolframalpha.com/input/?i=solve+for+t+d%2Fdt%5Bsqrt%28%28-z_1+t+%2B+t+v_1+%2B+x_1+-+y_1%29%5E2+%2B+%28-z_2+t+%2B+t+v_2+%2B+x_2+-+y_2%29%5E2+%2B+%28t+v_3+-+t+z_3+%2B+x_3+-+y_3%29%5E2%29%5D

  const auto &vi = i.getVelocity();
  const auto &vj = j.getVelocity();

  // Get old time step position
  const auto old_r_i = sub(i.getR(), mulScalar(vi, dt));
  const auto old_r_j = sub(j.getR(), mulScalar(vj, dt));

  // Compute nominator dot products
  const auto vi_ri = dot(vi, old_r_i);
  const auto vi_rj = dot(vi, old_r_j);
  const auto vj_ri = dot(vj, old_r_i);
  const auto vj_rj = dot(vj, old_r_j);

  const auto nominator = vi_rj + vj_ri - vi_ri - vj_rj;

  // Compute denominator dot products
  const auto two_vi_vj = 2.0 * dot(vi, vj);
  const auto vi_square = dot(vi, vi);
  const auto vj_square = dot(vj, vj);

  const auto denominator = vi_square + vj_square - two_vi_vj;

  // Compute t at minimal distance
  auto t = nominator / denominator;

  // If in the past, minimum is at t = 0
  // Else If in future timesteps, minimum for this is at t = dt
  t = std::clamp(t, 0., dt);

  // Compute actual distance by propagating along the line to t
  const auto p1 = add(old_r_i, mulScalar(vi, t));
  const auto p2 = add(old_r_j, mulScalar(vj, t));

  return std::make_tuple(sub(p1, p2), p2);
}
}  // namespace LADDS