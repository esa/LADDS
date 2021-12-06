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
#include "ladds/io/ConjunctionLogger.h"
#include "ladds/io/Timers.h"
#include "ladds/particle/Constellation.h"
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
   * @return Tuple of csvWriter, AccelerationAccumulator, and Integrator
   */
  [[nodiscard]] std::tuple<std::unique_ptr<FileOutput<AutoPas_t>>,
                           std::unique_ptr<Acceleration::AccelerationAccumulator<AutoPas_t>>,
                           std::unique_ptr<Integrator<AutoPas_t>>>
  initIntegrator(AutoPas_t &autopas, const YAML::Node &config);

  /**
   * Load the particles from the input csv files in the config as particles into AutoPas.
   *
   * @note Paths for csv files are relative to ladds/data!
   *
   * @param autopas
   * @param config
   */
  void loadSatellites(AutoPas_t &autopas, const YAML::Node &config);

  /**
   * Parse constellation information and prepare satellites for insertion.
   * @param config
   * @return Vector of Constellations
   */
  std::vector<Constellation> loadConstellations(const YAML::Node &config);

  /**
   * Tick constellation state machines and if applicable insert new satellites as well as delayed ones from previous
   * launch phase.
   * @param autopas
   * @param constellations
   * @param constellationInsertionFrequency
   * @param delayedInsertion Vector of satellites that could not be inserted in the last phase. This is an in/out
   * parameter!
   */
  void updateConstellation(AutoPas_t &autopas,
                           std::vector<Constellation> constellations,
                           int constellationInsertionFrequency,
                           std::vector<Particle> &delayedInsertion);

  /**
   * Check for collisions / conjunctions and write statistics about them.
   * @param autopas
   */
  void collisionDetection(size_t iteration,
                          Simulation::AutoPas_t &autopas,
                          ConjunctionLogger &conjunctionLogger,
                          size_t &totalConjunctions,
                          size_t progressOutputFrequency);

  /**
   * Write a vtk file with the current state of the simulation.
   */
  void writeVTK(size_t iteration, const AutoPas_t &autopas);

  /**
   * The main loop.
   * @param autopas
   * @param integrator
   * @param config
   */
  void simulationLoop(AutoPas_t &autopas,
                      Integrator<AutoPas_t> &integrator,
                      std::vector<Constellation> &constellations,
                      const YAML::Node &config);

  /**
   * Inserts particles into autopas if they have a safe distance to existing particles.
   * @param autopas Container where particles are about to be added.
   * @param newSatellites Satellites to be added.
   * @param cutoff
   * @return Vector of particles that were not added.
   */
  std::vector<Particle> checkedInsert(autopas::AutoPas<Particle> &autopas,
                                      const std::vector<Particle> &newSatellites,
                                      double cutoff);

  /**
   * One logger to log them all.
   */
  Logger &logger;

  /**
   * All timers used throughout the simulation.
   */
  Timers timers{};
};
