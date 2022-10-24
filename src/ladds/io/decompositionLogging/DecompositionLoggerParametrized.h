/**
 * @file DecompositionLoggerParametrized.h
 * @author F. Gratl
 * @date 24.10.22
 */

#pragma once

#include "DecompositionLogger.h"

namespace LADDS {

/**
 * This class serves as an intermediate helper interface to allow child classes to make use of CRTP while still allowing
 * a template free main interface via DecompositionLogger.
 * @tparam DecompositionT
 */
template <class DecompositionT>
class DecompositionLoggerParametrized : public DecompositionLogger {
 public:
  DecompositionLoggerParametrized(const std::string &loggerName,
                                  ConfigReader &config,
                                  std::string fileExtension,
                                  const DecompositionT &decomposition);

 protected:
  /**
   * Reference to the domain decomposition.
   */
  const DecompositionT &decomposition;

  /**
   * Generates a filename for the payload file (vts) of a specific rank.
   * @param iteration
   * @param rank If omitted, the current rank number is chosen.
   * @return
   */
  std::string filenamePayload(size_t iteration, int rank = -1) const;

  /**
   * Generates a filename for the metadata file (pvts).
   * @param iteration
   * @return
   */
  std::string filenameMetadata(size_t iteration) const;
};

template <class DecompositionT>
DecompositionLoggerParametrized<DecompositionT>::DecompositionLoggerParametrized(const std::string &loggerName,
                                                                                 ConfigReader &config,
                                                                                 std::string fileExtension,
                                                                                 const DecompositionT &decomposition)
    : DecompositionLogger(loggerName, config, fileExtension), decomposition(decomposition) {}

template <class DecompositionT>
std::string DecompositionLoggerParametrized<DecompositionT>::filenamePayload(size_t iteration, int rank) const {
  if (rank == -1) {
    autopas::AutoPas_MPI_Comm_rank(decomposition.getCommunicator(), &rank);
  }
  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(decomposition.getCommunicator(), &numRanks);
  std::stringstream ss;
  ss << "Decomp_" << std::setfill('0') << std::setw(static_cast<int>(std::to_string(numRanks).size())) << rank << "_"
     << std::setw(static_cast<int>(maxDigitsIterations)) << iteration << "." << this->fileExtension;
  return ss.str();
}

template <class DecompositionT>
std::string LADDS::DecompositionLoggerParametrized<DecompositionT>::filenameMetadata(size_t iteration) const {
  std::stringstream ss;
  ss << "Decomp_" << std::setfill('0') << std::setw(static_cast<int>(maxDigitsIterations)) << iteration << ".p"
     << this->fileExtension;
  return ss.str();
}

}  // namespace LADDS