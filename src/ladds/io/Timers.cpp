/**
 * @file Timers.cpp
 * @author F. Gratl
 * @date 06.12.21
 */

#include "Timers.h"

#include <iostream>

#include "ConfigReader.h"
#include "ladds/distributedMemParallelization/DomainDecomposition.h"

namespace LADDS {

void Timers::printTimers(ConfigReader &config, const DomainDecomposition &decomp) const {
  // local helper function to sum all mpi timers
  auto accMPITimers = [&](const autopas::utils::Timer &timer) {
    const auto localTime = timer.getTotalTime();
    long reducedTime{};
    autopas::AutoPas_MPI_Reduce(
        &localTime, &reducedTime, 1, AUTOPAS_MPI_LONG, AUTOPAS_MPI_SUM, 0l, decomp.getCommunicator());
    return reducedTime;
  };

  const auto iterations = config.get<size_t>("sim/iterations");
  const auto timeTotalAcc = accMPITimers(total);
  const auto timeSimAcc = accMPITimers(simulation);
  const auto timeTotal = total.getTotalTime();
  const auto timeSim = simulation.getTotalTime();
  const auto maxNumberOfDigits = static_cast<int>(std::to_string(timeTotalAcc).length());

  // print imbalance
  auto printLB = [&](const auto &timer) {
    // only print timers that were actually used. This assumes that no timer can run on only one rank.
    if (timer.getTotalTime() == 0) {
      return std::string();
    }

    auto [avgImbalance, maxImbalance] = calcImbalances(timer, decomp);

    std::stringstream ss;
    ss << std::setprecision(floatStringPrecision);
    ss << " | relative LB mean: " << std::setw(floatStringPrecision + 6) << avgImbalance << " max: " << maxImbalance
       << "\n";
    return ss.str();
  };

  SPDLOG_LOGGER_INFO(
      config.getLogger().get(),
      "Runtime summary:\n{}",
      // clang-format off
      timerToString("Total (ranks accumulated)          ", timeTotalAcc, maxNumberOfDigits) + printLB(total) +
      timerToString("  Initialization                   ", accMPITimers(initialization), maxNumberOfDigits, timeTotalAcc) +
          printLB(initialization) +
      timerToString("  Simulation                       ", timeSimAcc, maxNumberOfDigits, timeTotalAcc) +
          printLB(simulation) +
      timerToString("    Integrator                     ", accMPITimers(integrator), maxNumberOfDigits, timeSimAcc) +
          printLB(integrator) +
      timerToString("    Resolving Burn ups             ", accMPITimers(burnUps), maxNumberOfDigits, timeSimAcc) +
          printLB(burnUps) +
      timerToString("    Constellation insertion        ", accMPITimers(constellationInsertion), maxNumberOfDigits, timeSimAcc) +
          printLB(constellationInsertion) +
      timerToString("    Communication                  ", accMPITimers(particleCommunication), maxNumberOfDigits, timeSimAcc) +
          printLB(particleCommunication) +
      timerToString("    Collision detection            ", accMPITimers(collisionDetection), maxNumberOfDigits, timeSimAcc) +
          printLB(collisionDetection) +
      timerToString("    Collision detection immigrants ", accMPITimers(collisionDetectionImmigrants), maxNumberOfDigits, timeSimAcc) +
          printLB(collisionDetectionImmigrants) +
      timerToString("    Collision detection emmigrants ", accMPITimers(collisionDetectionEmmigrants), maxNumberOfDigits, timeSimAcc) +
          printLB(collisionDetectionEmmigrants) +   
      timerToString("    Collision writer               ", accMPITimers(collisionWriting), maxNumberOfDigits, timeSimAcc) +
          printLB(collisionWriting) +
      timerToString("    Collision simulation           ", accMPITimers(collisionSimulation), maxNumberOfDigits, timeSimAcc) +
          printLB(collisionSimulation) +
      timerToString("    Collision writer               ", accMPITimers(evasionWriting), maxNumberOfDigits, timeSimAcc) +
          printLB(evasionWriting) +
      timerToString("    Container update               ", accMPITimers(containerUpdate), maxNumberOfDigits, timeSimAcc) +
          printLB(containerUpdate) +
      timerToString("    Output                         ", accMPITimers(output), maxNumberOfDigits, timeTotalAcc) +
          printLB(output /*offest for formatting*/) +
      timerToString("Total (wall-time)                  ", timeTotal, maxNumberOfDigits) + "\n" +
      timerToString("  One iteration                    ", timeSim / iterations, maxNumberOfDigits, timeTotal));
  // clang-format on
}

/**
 * Turns the timers into a human readable string.
 * @param name: The timer's name.
 * @param timeNS: The time in nano seconds.
 * @param numberWidth: The precision of the printed number.
 * @param maxTime: The simulation's total execution time.
 * @return All information of the timer in a human readable string.
 *
 * @note Taken from md-flexible.
 */
std::string Timers::timerToString(const std::string &name, long timeNS, int numberWidth, long maxTime) {
  // only print timers that were actually used
  if (timeNS == 0) {
    return "";
  }

  std::ostringstream ss;
  ss << std::fixed << std::setprecision(floatStringPrecision) << name << " : " << std::setw(numberWidth) << std::right
     << timeNS
     << " ns ("
     // min width of the representation of seconds is numberWidth - 9 (from conversion) + 4 (for dot and digits after)
     << std::setw(numberWidth - 5) << ((double)timeNS * 1e-9) << "s)";
  if (maxTime != 0) {
    ss << " =" << std::setw(7) << std::right << ((double)timeNS / (double)maxTime * 100) << "%";
  } else {
    // fill up same space with whitespaces
    ss << std::setw(2 + 7 + 1) << "";
  }
  return ss.str();
}

std::tuple<double, double> Timers::calcImbalances(const autopas::utils::Timer &timer,
                                                  const DomainDecomposition &decomp) {
  const auto localTime = timer.getTotalTime();
  long sumOfTimes{0l};
  autopas::AutoPas_MPI_Allreduce(
      &localTime, &sumOfTimes, 1, AUTOPAS_MPI_LONG, AUTOPAS_MPI_SUM, decomp.getCommunicator());
  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(decomp.getCommunicator(), &numRanks);
  const auto avgTime = static_cast<double>(sumOfTimes) / numRanks;
  const auto localImbalance = (static_cast<double>(localTime) - avgTime) / avgTime;
  double sumImbalances{0.};
  autopas::AutoPas_MPI_Reduce(
      &localImbalance, &sumImbalances, 1, AUTOPAS_MPI_DOUBLE, AUTOPAS_MPI_SUM, 0l, decomp.getCommunicator());
  const double avgImbalance = sumImbalances / numRanks;
  double maxImbalance{0.};
  autopas::AutoPas_MPI_Reduce(
      &localImbalance, &maxImbalance, 1, AUTOPAS_MPI_DOUBLE, AUTOPAS_MPI_MAX, 0l, decomp.getCommunicator());
  return std::make_tuple(avgImbalance, maxImbalance);
}

}  // namespace LADDS