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

#include <memory>

#include "CollisionFunctor.h"
#include "TypeDefinitions.h"
#include "ladds/io/ConfigReader.h"
#include "ladds/io/ConjunctionLogger.h"
#include "ladds/io/Timers.h"
#include "ladds/io/hdf5/HDF5Writer.h"
#include "ladds/particle/Constellation.h"
#include "ladds/particle/Particle.h"

extern template class autopas::AutoPas<Particle>;
extern template bool autopas::AutoPas<Particle>::iteratePairwise(CollisionFunctor *);

/**
 * The main simulation class responsible for all high-level logic.
 */
class Simulation {
 public:
  /**
   * Constructor.
   * @param logger The class stores a reference to a logger but does not take ownership.
   */
  explicit Simulation(Logger &logger) : logger(logger) {}

  /**
   * Run the whole simulation with the given config.
   * @param config
   */
  void run(ConfigReader &config);

  /**
   * Create and initialize an AutoPas object from the given config.
   * @param config
   * @return
   */
  [[nodiscard]] std::unique_ptr<AutoPas_t> initAutoPas(ConfigReader &config);

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
  initIntegrator(AutoPas_t &autopas, ConfigReader &config);

  /**
   * Tick constellation state machines and if applicable insert new satellites as well as delayed ones from previous
   * launch phase.
   * @param autopas
   * @param constellations
   * @param delayedInsertion Vector of satellites that could not be inserted in the last phase. This is an in/out
   * parameter!
   * @param constellationCutoff range parameter for checked insertion: if the insertion would be within a distance
   * of constellationCutoff to any other object the insertion is delayed instead
   */
  void updateConstellation(AutoPas_t &autopas,
                           std::vector<Constellation> &constellations,
                           std::vector<Particle> &delayedInsertion,
                           double constellationCutoff);

  /**
   * Check for collisions / conjunctions and write statistics about them.
   * @param autopas
   */
  std::unordered_map<Particle *, std::tuple<Particle *, double>> collisionDetection(AutoPas_t &autopas,
                                                                                    double deltaT,
                                                                                    double conjunctionThreshold);

  /**
   * The main loop.
   * @param autopas
   * @param integrator
   * @param config
   */
  void simulationLoop(AutoPas_t &autopas,
                      Integrator<AutoPas_t> &integrator,
                      std::vector<Constellation> &constellations,
                      ConfigReader &config);

  /**
   * Inserts particles into autopas if they have a safe distance to existing particles.
   * @param autopas Container where particles are about to be added.
   * @param newSatellites Satellites to be added.
   * @param constellationCutoff sphere with radius constellationCutoff around inserted satellite must be
   * empty for satellite to be inserted or else the insertion will be delayed (passed in the returned vector)
   * @return Vector of particles that were not added.
   */
  std::vector<Particle> checkedInsert(autopas::AutoPas<Particle> &autopas,
                                      const std::vector<Particle> &newSatellites,
                                      double constellationCutoff);

  /**
   * One logger to log them all.
   */
  Logger &logger;

  /**
   * All timers used throughout the simulation.
   */
  Timers timers{};
};
