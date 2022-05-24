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
#include <satellitePropagator/physics/YoshidaIntegrator.h>

#include <memory>

#include "CollisionFunctor.h"
#include "TypeDefinitions.h"
#include "ladds/io/ConfigReader.h"
#include "ladds/io/ConjunctionLogger.h"
#include "ladds/io/Timers.h"
#include "ladds/io/hdf5/HDF5Writer.h"
#include "ladds/particle/Constellation.h"
#include "ladds/particle/Particle.h"
#include "ladds/distributedMemParallelization/DomainDecomposition.h"

namespace LADDS {

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
  [[nodiscard]] std::unique_ptr<AutoPas_t> initAutoPas(ConfigReader &config, DomainDecomposition &domainDecomp);

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
                           std::unique_ptr<YoshidaIntegrator<AutoPas_t>>>
  initIntegrator(AutoPas_t &autopas, ConfigReader &config);

  /**
   * Depending on config initialize readers.
   * @param config
   * @return tuple<hdf5WriteFrequency, hdf5Writer, conjuctionWriter>
   */
  [[nodiscard]] std::tuple<size_t, std::shared_ptr<HDF5Writer>, std::shared_ptr<ConjuctionWriterInterface>> initWriter(
      ConfigReader &config);

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
   * @param deltaT
   * @param collisionDistanceFactor See CollisionFunctor::_collisionDistanceFactor
   * @param minDetectionRadius
   * @return Tuple of the collisions and whether AutoPas is currently in tuning mode.
   */
  std::tuple<CollisionFunctor::CollisionCollectionT, bool> collisionDetection(AutoPas_t &autopas,
                                                                              double deltaT,
                                                                              double collisionDistanceFactor,
                                                                              double minDetectionRadius);

  /**
   * Updates the configuration with the latest AutoPas configuration and writes it to a new YAML file.
   * @param config
   * @param autopas
   */
  void dumpCalibratedConfig(ConfigReader &config, const AutoPas_t &autopas) const;

  /**
   * The main loop.
   * @param autopas
   * @param integrator
   * @param constellations
   * @param config
   * @param domainDecomposition
   * @return Number of observed collisions.
   */
  [[maybe_unused]] size_t simulationLoop(AutoPas_t &autopas,
                                         YoshidaIntegrator<AutoPas_t> &integrator,
                                         std::vector<Constellation> &constellations,
                                         ConfigReader &config,
                                         DomainDecomposition &domainDecomposition);

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
   * Distributes global IDs for all constellation particles and returns the resulting maxExistingParticleId.
   * @param autopas
   * @param constellations
   * @return size_t maxExistingParticleId
   */
  size_t setConstellationIDs(autopas::AutoPas<Particle> &autopas, std::vector<Constellation> &constellations);

  /**
   * Remove all particles below a certain altitude from the particle container.
   * @param autopas
   * @param burnUpAltitude Height above ground. [km]
   */
  void deleteBurnUps(autopas::AutoPas<Particle> &autopas, double burnUpAltitude) const;

  /**
   * Computes a timeout value in seconds from the information given in the config. If nothing is given in the config
   * the function returns 0.
   * @note A timeout of 0 is considered to be no timeout
   * @param config
   * @return timeout in seconds.
   */
  size_t computeTimeout(ConfigReader &config);
  /**
   * One logger to log them all.
   */
  Logger &logger;

  /**
   * All timers used throughout the simulation.
   */
  Timers timers{};
};

}  // namespace LADDS
