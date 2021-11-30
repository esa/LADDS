/**
 * @file Simulation.h
 * @author F. Gratl
 * @date 30.11.21
 */

#pragma once

#include <autopas/AutoPasDecl.h>
#include <satellitePropagator/io/FileOutput.h>
#include <satellitePropagator/physics/AccelerationAccumulator.h>
#include <satellitePropagator/physics/Integrator.h>
#include <yaml-cpp/yaml.h>

#include "CollisionFunctor.h"
#include "ladds/particle/Particle.h"

extern template class autopas::AutoPas<Particle>;
extern template bool autopas::AutoPas<Particle>::iteratePairwise(CollisionFunctor *);

class Simulation {
 public:
  using AutoPas_t = autopas::AutoPas<Particle>;

  void run(const YAML::Node &config);

 private:
  [[nodiscard]] std::unique_ptr<AutoPas_t> initAutoPas(const YAML::Node &config);
  [[nodiscard]] std::unique_ptr<Integrator<Simulation::AutoPas_t>> initIntegrator(AutoPas_t &autopas,
                                                                                  const YAML::Node &config);
  void loadSatellites(AutoPas_t &autopas, const YAML::Node &config);

  void simulationLoop(AutoPas_t &autopas, Integrator<AutoPas_t> &integrator, const YAML::Node &config);

  void printTimers(const YAML::Node &config);

  std::string timerToString(const std::string &name, long timeNS, int numberWidth = 0l, long maxTime = 0ul);

  /**
   * Floating point precision for command line output.
   */
  static constexpr int floatStringPrecision = 3;

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
