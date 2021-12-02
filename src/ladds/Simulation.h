/**
 * @file Simulation.h
 * @author F. Gratl
 * @date 30.11.21
 */

#pragma once

#include <autopas/AutoPasDecl.h>
#include <ladds/io/Logger.h>
#include <satellitePropagator/io/FileOutput.h>
#include <satellitePropagator/physics/AccelerationAccumulator.h>
#include <satellitePropagator/physics/Integrator.h>
#include <yaml-cpp/yaml.h>

#include <memory>

#include "CollisionFunctor.h"
#include "ladds/particle/Particle.h"

extern template class autopas::AutoPas<Particle>;
extern template bool autopas::AutoPas<Particle>::iteratePairwise(CollisionFunctor *);

/**
 * The main simulation class responsible for all high-level logic.
 */
class Simulation {
 public:
  using AutoPas_t = autopas::AutoPas<Particle>;

  /**
   * Constructor.
   * @param logger The class stores a reference to a logger but does not take ownership.
   */
  explicit Simulation(Logger &logger) : logger(logger) {}

  /**
   * Run the whole simulation with the given config.
   * @param config
   */
  void run(const YAML::Node &config);

 private:
  /**
   * Create and initialize an AutoPas object from the given config.
   * @param config
   * @return
   */
  [[nodiscard]] std::unique_ptr<AutoPas_t> initAutoPas(const YAML::Node &config);

  /**
   * Create and initialize the integrator.
   * As the integrator internally needs the csvWriter and AccelerationAccumulator but does not take ownership
   * these objects. Hence, they need to be kept alive as long as the integrator is needed.
   * @param autopas
   * @param config
   * @return
   */
  [[nodiscard]] auto initIntegrator(AutoPas_t &autopas, const YAML::Node &config);

  /**
   * Load the particles from the input csv files in the config as particles into AutoPas.
   * @param autopas
   * @param config
   */
  void loadSatellites(AutoPas_t &autopas, const YAML::Node &config);

  /**
   * The main loop.
   * @param autopas
   * @param integrator
   * @param config
   */
  void simulationLoop(AutoPas_t &autopas, Integrator<AutoPas_t> &integrator, const YAML::Node &config);

  /**
   * Print timer information to stdout.
   * @param config
   */
  void printTimers(const YAML::Node &config) const;

  /**
   * Helper function to pretty-print timers.
   * @param name
   * @param timeNS
   * @param numberWidth
   * @param maxTime
   * @return
   */
  static std::string timerToString(const std::string &name, long timeNS, int numberWidth = 0, long maxTime = 0ul);

  /**
   * Floating point precision for command line output.
   */
  static constexpr int floatStringPrecision = 3;

  /**
   * One logger to log them all.
   */
  Logger &logger;

  /**
   * All timers used throughout the simulation.
   */
  struct Timers {
    autopas::utils::Timer total{};
    autopas::utils::Timer initialization{};
    autopas::utils::Timer simulation{};
    autopas::utils::Timer integrator{};
    autopas::utils::Timer collisionDetection{};
    autopas::utils::Timer containerUpdate{};
    autopas::utils::Timer output{};
  } __attribute__((aligned(128))) timers{};
};
