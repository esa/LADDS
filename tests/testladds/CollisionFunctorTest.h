/**
 * @file CollisionFunctorTest.h
 * @author F. Gratl
 * @date 05/07/2021
 */

#pragma once

#include <gtest/gtest.h>
#include <autopas/utils/ArrayUtils.h>

using ParameterTuple =
    std::tuple<std::array<double, 3>,std::array<double, 3>,std::array<double, 3>,std::array<double, 3>,double,double>;

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
      const auto& [x1,x2,v1,v2,dt,expected_dist] = static_cast<ParamType>(info.param);

      auto inputTuple = static_cast<ParamType>(info.param);
      std::stringstream ss;
      ss << "x1_" << autopas::utils::ArrayUtils::to_string(x1, "", {"",""}) << "x2_" << autopas::utils::ArrayUtils::to_string(x2, "", {"",""}) << "_"
         << "v1_" << autopas::utils::ArrayUtils::to_string(v1, "", {"",""}) << "v2_" << autopas::utils::ArrayUtils::to_string(v2, "", {"",""}) << "_"
         << "dt_" << dt << "_" << "expected_dist_" << expected_dist;
      auto str = ss.str();
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      std::replace(str.begin(), str.end(), ' ', '_');
      return str;
    }
  };
};