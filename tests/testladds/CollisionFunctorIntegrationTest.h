/**
 * @file CollisionFunctorTest.h
 * @author F. Gratl
 * @date 05/07/2021
 */

#pragma once
#include <autopas/options/ContainerOption.h>
#include <autopas/options/DataLayoutOption.h>
#include <autopas/options/Newton3Option.h>
#include <autopas/options/TraversalOption.h>
#include <gtest/gtest.h>
#include <ladds/Debris.h>

#include "tuple"
#include "vector"

using ParameterTuple =
    std::tuple<autopas::TraversalOption, autopas::DataLayoutOption, autopas::Newton3Option, double /*cellSizeFactor*/>;

class CollisionFunctorIntegrationTest : public testing::TestWithParam<ParameterTuple> {
 public:
  static void SetUpTestSuite();

  /**
   * Helper struct for pretty printing the generated test names.
   */
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      const auto &[traversal, dataLayout, newton3, cellSizeFactor] = static_cast<ParamType>(info.param);

      auto inputTuple = static_cast<ParamType>(info.param);
      std::stringstream ss;
      ss << traversal << "_" << dataLayout << "_"
         << "N3_" << newton3 << "_"
         << "CSF_" << cellSizeFactor;
      auto str = ss.str();
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      return str;
    }
  };

  static inline std::vector<std::pair<size_t, size_t>> _reference{};
  static inline std::vector<Debris> _debris{};
  static constexpr double _cutoff{1.};
  static constexpr std::array<double, 3> _boxMin{0., 0., 0.};
  static constexpr std::array<double, 3> _boxMax{5., 5., 5.};
};
