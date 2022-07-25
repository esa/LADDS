/**
 * @file CollisionFunctorTest.h
 * @author F. Gratl
 * @date 05/07/2021
 */

#pragma once

#include <autopas/utils/ArrayUtils.h>
#include <gtest/gtest.h>

using ParameterTuple = std::tuple<std::array<double, 3>,
                                  std::array<double, 3>,
                                  std::array<double, 3>,
                                  std::array<double, 3>,
                                  double,
                                  double,
                                  std::array<double, 3>>;

/**
 * Simple test for basic functionalities of the CollisionFunctor.
 */
class CollisionFunctorTest : public testing::TestWithParam<ParameterTuple> {
 public:
  /**
   * Helper struct for pretty printing the generated test names.
   */
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      const auto &[x1, x2, v1, v2, dt, expected_dist, collisionPoint] = static_cast<ParamType>(info.param);

      std::stringstream ss;
      ss << "x1" << autopas::utils::ArrayUtils::to_string(x1, "", {"_", "_"}) << "x2"
         << autopas::utils::ArrayUtils::to_string(x2, "", {"_", "_"}) << "v1"
         << autopas::utils::ArrayUtils::to_string(v1, "", {"_", "_"}) << "v2"
         << autopas::utils::ArrayUtils::to_string(v2, "", {"_", "_"}) << "dt_" << dt << "_"
         << "expected_dist_" << expected_dist << "collisionPoint"
         << autopas::utils::ArrayUtils::to_string(collisionPoint, "", {"_", "_"});
      auto str = ss.str();
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      std::replace(str.begin(), str.end(), ' ', '_');
      return str;
    }
  };
};