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
#include "ladds/distributedMemParallelization/DomainDecomposition.h"
#include "ladds/distributedMemParallelization/RegularGridDecomposition.h"
#include "ladds/io/ConfigReader.h"
#include "ladds/io/ConjunctionLogger.h"
#include "ladds/io/Timers.h"
#include "ladds/io/hdf5/HDF5Writer.h"
#include "ladds/particle/BreakupWrapper.h"
#include "ladds/particle/Constellation.h"
#include "ladds/particle/Particle.h"

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
   * Interact all incoming particles with all particles which potentially crossed its path since the last container
   * update.
   *
   * These particles are found in a box around the immigrant's position -deltaT time ago.
   * The box has a side length of 2x the maximum coverable distance by any particle.
   *
   * @note The first pointers in the returned tuple collection point to particles in the immigrant vector!
   *
   * @param autopas
   * @param incomingParticles
   * @param deltaT Time since the last container update.
   * @param maxV Maximal velocity a particle is assumed to have. Has to be positive.
   * @param collisionDistanceFactor See CollisionFunctor::_collisionDistanceFactor
   * @param minDetectionRadius
   * @return Collection of collision partners
   */
  CollisionFunctor::CollisionCollectionT collisionDetectionImmigrants(AutoPas_t &autopas,
                                                                      std::vector<Particle> &incomingParticles,
                                                                      double deltaT,
                                                                      double maxV,
                                                                      double collisionDistanceFactor,
                                                                      double minDetectionRadius);

  /**
   * Auxiliary function to avoid code duplication. Triggers the writing of collisions and if necessary the breakup
   * simulation. Also invokes the relevant timers.
   * @param iteration
   * @param collisions
   * @param conjuctionWriter
   * @param breakupWrapper
   */
  void processCollisions(size_t iteration,
                         const CollisionFunctor::CollisionCollectionT &collisions,
                         ConjuctionWriterInterface &conjunctionWriter,
                         BreakupWrapper *breakupWrapper);

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
   * @return
   */
  void setConstellationIDs(autopas::AutoPas<Particle> &autopas, std::vector<Constellation> &constellations);

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
   * Send the given list of leaving particles to all (up to) 26 logical surrounding ranks and receive their leaving
   * particles which are relevant for the local rank.
   * @param leavingParticles in/out parameter of leaving particles. If everything worked the vector should be empty
   * after the function call.
   * @param autopas
   * @param decomposition
   * @return Vector of incoming particles.
   */
  std::vector<Particle> communicateParticles(std::vector<Particle> &leavingParticles,
                                             autopas::AutoPas<Particle> &autopas,
                                             const RegularGridDecomposition &decomposition);

  /**
   * Prints one line on the info log level stating every rank's number of particles, sorted by rank number.
   * @param autopas
   * @param decomposition
   */
  void printNumParticlesPerRank(const autopas::AutoPas<Particle> &autopas,
                                const DomainDecomposition &decomposition) const;

  /**
   * Loggs global information about simulation progress to std::out from rank 1.
   * @note Contains global communication to obtain full information.
   * @param iteration
   * @param numParticlesLocal
   * @param totalConjunctionsLocal
   * @param comm
   */
  void printProgressOutput(size_t iteration,
                           size_t numParticlesLocal,
                           size_t totalConjunctionsLocal,
                           const autopas::AutoPas_MPI_Comm &comm);
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
