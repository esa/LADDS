/**
 * @file DecompositionLogger.h
 * @author F. Gratl
 * @date 24.05.22
 */

#pragma once

namespace LADDS {
class DecompositionLogger {
 public:
  virtual ~DecompositionLogger() = default;
  /**
   * Creates the pvts file that references all vts files
   * @note Should only be called by one rank!
   * @param iteration
   */
  virtual void writeMetafile(size_t iteration) const = 0;

  /**
   * Creates all vts files which contain the actual decomposition.
   * @param iteration
   * @param autoPasConfig
   */
  virtual void writePayload(size_t iteration, const autopas::Configuration &autoPasConfig) const = 0;
};
}  // namespace LADDS