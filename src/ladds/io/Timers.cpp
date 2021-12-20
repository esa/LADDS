/**
 * @file Timers.cpp
 * @author F. Gratl
 * @date 06.12.21
 */

#include "Timers.h"

#include <iostream>

#include "ConfigReader.h"

void Timers::printTimers(ConfigReader &config) const {
  const auto iterations = config.get<size_t>("sim/iterations");
  const auto timeTotal = total.getTotalTime();
  const auto timeSim = simulation.getTotalTime();
  const auto maximumNumberOfDigits = static_cast<int>(std::to_string(timeTotal).length());
  std::cout << timerToString("Total                       ", timeTotal, maximumNumberOfDigits);
  std::cout << timerToString(
      "  Initialization            ", initialization.getTotalTime(), maximumNumberOfDigits, timeTotal);
  std::cout << timerToString("  Simulation                ", timeSim, maximumNumberOfDigits, timeTotal);
  std::cout << timerToString("    Integrator              ", integrator.getTotalTime(), maximumNumberOfDigits, timeSim);
  std::cout << timerToString(
      "    Constellation insertion ", constellationInsertion.getTotalTime(), maximumNumberOfDigits, timeSim);
  std::cout << timerToString(
      "    Collision detection     ", collisionDetection.getTotalTime(), maximumNumberOfDigits, timeSim);
  std::cout << timerToString(
      "    Container update        ", containerUpdate.getTotalTime(), maximumNumberOfDigits, timeSim);
  std::cout << timerToString("    Output                  ", output.getTotalTime(), maximumNumberOfDigits, timeTotal);
  std::cout << timerToString("One iteration               ", timeSim / iterations, maximumNumberOfDigits, timeTotal);
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
  }
  ss << std::endl;
  return ss.str();
}