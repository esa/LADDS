/**
 * @file Simulation.cpp
 * @author F. Gratl
 * @date 30.11.21
 */
#include "Simulation.h"

#include <breakupModel/output/VTKWriter.h>

// needed for autopas::AutoPas::getCurrentConfiguration()
#include <autopas/AutoPasImpl.h>

#include <array>
#include <limits>
#include <tuple>
#include <vector>

#include "ladds/distributedMemParallelization/AltitudeBasedDecomposition.h"
#include "ladds/distributedMemParallelization/RankMigration.h"
#include "ladds/io/ConjunctionLogger.h"
#include "ladds/io/SatelliteLoader.h"
#include "ladds/io/VTUWriter.h"
#include "ladds/io/decompositionLogging/DecompositionLogger.h"
#include "ladds/io/decompositionLogging/RegularGridDecompositionLogger.h"
#include "ladds/io/hdf5/HDF5Writer.h"
#include "ladds/particle/Constellation.h"

// Declare the main AutoPas class as extern template instantiation. It is instantiated in AutoPasClass.cpp.
extern template class autopas::AutoPas<LADDS::Particle>;
extern template bool autopas::AutoPas<LADDS::Particle>::iteratePairwise(LADDS::CollisionFunctor *);

namespace LADDS {

/**
 * Anonymous namespace to hide this content from the outside.
 */
namespace {
/**
 * Local helper function for parsing AutoPas Set options or setting default values.
 * @tparam Option
 * @tparam F
 * @param node
 * @param setterFun
 * @param defaultVals
 */
template <class Option, class F>
void setAutoPasOption(ConfigReader &config,
                      const std::string &valuePath,
                      F setterFun,
                      const std::set<Option> &defaultVals) {
  auto valStr = config.template get<std::string>(
      valuePath, autopas::utils::ArrayUtils::to_string(defaultVals, " ", {"", ""}), true);
  const auto options = Option::parseOptions(valStr);
  setterFun(options);
}
}  // namespace

std::unique_ptr<AutoPas_t> Simulation::initAutoPas(ConfigReader &config, DomainDecomposition &domainDecomp) {
  auto autopas = std::make_unique<AutoPas_t>();

  const auto maxAltitude = config.get<double>("sim/maxAltitude");
  const auto cutoff = config.get<double>("autopas/cutoff");
  const auto desiredCellsPerDimension = config.get<double>("autopas/desiredCellsPerDimension", 25);
  const auto deltaT = config.get<double>("sim/deltaT");
  const auto verletRebuildFrequency = config.get<unsigned int>("autopas/rebuildFrequency", 1);
  // *8.5 : Assumed max km/s   TODO: get this from input?
  // *2: because particles might flight directly towards each other
  const auto verletSkin = 2 * 8.5 * deltaT * verletRebuildFrequency;

  SPDLOG_LOGGER_DEBUG(logger.get(), "Verlet Skin: {}", verletSkin);

  SPDLOG_LOGGER_DEBUG(
      logger.get(), "Setting BoxMin / Max: {} {}", domainDecomp.getLocalBoxMin(), domainDecomp.getLocalBoxMax());
  autopas->setBoxMin(domainDecomp.getLocalBoxMin());
  autopas->setBoxMax(domainDecomp.getLocalBoxMax());
  autopas->setCutoff(cutoff);
  autopas->setVerletSkin(verletSkin);
  autopas->setVerletRebuildFrequency(verletRebuildFrequency);
  // Scale Cell size so that we get the desired number of cells
  // -2 because internally there will be two halo cells added on top of maxAltitude
  // FIXME: adapt calculation of CSF to MPI
  autopas->setCellSizeFactor((maxAltitude * 2.) / ((cutoff + verletSkin) * (desiredCellsPerDimension - 2)));

  // hardcode values that seem to be optimal
  std::set<autopas::Newton3Option> optimalNewton3Opt{autopas::Newton3Option::enabled};
  std::set<autopas::DataLayoutOption> optimalDataLayoutOpt{autopas::DataLayoutOption::aos};
  std::set<autopas::ContainerOption> optimalContainerOpt{autopas::ContainerOption::linkedCells};
  std::set<autopas::TraversalOption> optimalTraversalOpt{autopas::TraversalOption::lc_c04_HCP};
  // if tuning mode is enabled use full range of AutoPas defaults
  if (config.get<bool>("autopas/tuningMode", false)) {
    optimalNewton3Opt = autopas->getAllowedNewton3Options();
    optimalDataLayoutOpt = autopas->getAllowedDataLayouts();
    optimalContainerOpt = autopas->getAllowedContainers();
    optimalTraversalOpt = autopas->getAllowedTraversals();
  }
  // however, always allow setting config options via yaml
  setAutoPasOption<autopas::Newton3Option>(
      config, "autopas/Newton3", [&](const auto &op) { autopas->setAllowedNewton3Options(op); }, optimalNewton3Opt);
  setAutoPasOption<autopas::DataLayoutOption>(
      config, "autopas/DataLayout", [&](const auto &op) { autopas->setAllowedDataLayouts(op); }, optimalDataLayoutOpt);
  setAutoPasOption<autopas::ContainerOption>(
      config, "autopas/Container", [&](const auto &op) { autopas->setAllowedContainers(op); }, optimalContainerOpt);
  setAutoPasOption<autopas::TraversalOption>(
      config, "autopas/Traversal", [&](const auto &op) { autopas->setAllowedTraversals(op); }, optimalTraversalOpt);
  // arbitrary number. Can be changed to whatever makes sense.
  autopas->setTuningInterval(std::numeric_limits<unsigned int>::max());
  autopas->setSelectorStrategy(autopas::SelectorStrategyOption::fastestMean);
  autopas->setNumSamples(verletRebuildFrequency);
  autopas::Logger::get()->set_level(spdlog::level::from_str(config.get<std::string>("autopas/logLevel", "error")));
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);
  autopas->setOutputSuffix("Rank" + std::to_string(rank) + "_");
  autopas->init();

  return autopas;
}

std::tuple<std::unique_ptr<FileOutput<AutoPas_t>>,
           std::unique_ptr<Acceleration::AccelerationAccumulator<AutoPas_t>>,
           std::unique_ptr<YoshidaIntegrator<AutoPas_t>>>
Simulation::initIntegrator(AutoPas_t &autopas, ConfigReader &config) {
  // initialization of the integrator
  std::array<bool, 8> selectedPropagatorComponents{};
  selectedPropagatorComponents = {config.get<bool>("sim/prop/useKEPComponent", false),
                                  config.get<bool>("sim/prop/useJ2Component", false),
                                  config.get<bool>("sim/prop/useC22Component", false),
                                  config.get<bool>("sim/prop/useS22Component", false),
                                  config.get<bool>("sim/prop/useSOLComponent", false),
                                  config.get<bool>("sim/prop/useLUNComponent", false),
                                  config.get<bool>("sim/prop/useSRPComponent", false),
                                  config.get<bool>("sim/prop/useDRAGComponent", false)};

  std::unique_ptr<FileOutput<AutoPas_t>> csvWriter{nullptr};
  if (const auto csvFilename = config.get<std::string>("io/csv/propagatorOutput", ""); not csvFilename.empty()) {
    csvWriter =
        std::make_unique<FileOutput<AutoPas_t>>(autopas, csvFilename, OutputFile::CSV, selectedPropagatorComponents);
  }
  auto accumulator = std::make_unique<Acceleration::AccelerationAccumulator<AutoPas_t>>(
      selectedPropagatorComponents, autopas, 0.0, csvWriter.get());
  auto deltaT = config.get<double>("sim/deltaT");
  auto integrator = std::make_unique<YoshidaIntegrator<AutoPas_t>>(autopas, *accumulator, deltaT);

  return std::make_tuple<>(std::move(csvWriter), std::move(accumulator), std::move(integrator));
}

std::tuple<size_t, std::shared_ptr<HDF5Writer>, std::shared_ptr<ConjuctionWriterInterface>> Simulation::initWriter(
    ConfigReader &config) {
  if (config.defines("io/hdf5")) {
    std::shared_ptr<HDF5Writer> hdf5Writer;
    const auto hdf5WriteFrequency = config.get<unsigned int>("io/hdf5/writeFrequency", 0u);
    const auto hdf5CompressionLvl = config.get<unsigned int>("io/hdf5/compressionLevel", 0u);
    // either append to checkpoint ...
    const auto checkpointFileName = config.get<std::string>("io/hdf5/checkpoint/file", "");
    if (not checkpointFileName.empty()) {
      // TODO: Remove DATADIR functionality
      const auto checkpointPath = std::string(DATADIR) + checkpointFileName;
      // compression level already set when file already exists
      hdf5Writer = std::make_shared<HDF5Writer>(checkpointPath, false, 0);
    }
    // ... or write to new HDF5 file ....
    if (const auto hdf5FileName =
            [&]() {
              int rank{};
              autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);
              const auto f = config.get<std::string>("io/hdf5/fileName", "");
              return f.empty() ? "" : f + "_rank_" + std::to_string(rank);
            }();
        not hdf5FileName.empty()) {
      // ... but not both!
      if (not checkpointFileName.empty()) {
        throw std::runtime_error(
            "HDF5Writer already defined!\nProbably both a checkpoint and a HDF5 output file are defined. Please choose "
            "only one.");
      }
      hdf5Writer = std::make_shared<HDF5Writer>(hdf5FileName, true, hdf5CompressionLvl);
    }
    // check that at least something set a filename
    if (not hdf5Writer) {
      throw std::runtime_error(
          "Config suggests HDF5 should be used (found io/hdf5) but neither a fileName nor a checkpoint to write to is "
          "defined.");
    }
    return {hdf5WriteFrequency, hdf5Writer, hdf5Writer};
  } else {
    return {0u, nullptr, std::make_shared<ConjunctionLogger>("")};
  }
}

void Simulation::updateConstellation(AutoPas_t &autopas,
                                     std::vector<Constellation> &constellations,
                                     std::vector<Particle> &delayedInsertionTotal,
                                     double constellationCutoff) {
  // first insert delayed particles from previous insertion and collect the repeatedly delayed
  delayedInsertionTotal = checkedInsert(autopas, delayedInsertionTotal, constellationCutoff);
  // container collecting delayed particles from one constellation at a time in order to append them to
  // totalDelayedInsertion
  std::vector<Particle> delayedInsertion;
  for (auto &constellation : constellations) {
    // new satellites are gradually added to the simulation according to their starting time and operation duration
    auto newSatellites = constellation.tick();
    delayedInsertion = checkedInsert(autopas, newSatellites, constellationCutoff);
    delayedInsertionTotal.insert(delayedInsertionTotal.end(), delayedInsertion.begin(), delayedInsertion.end());
  }
}

std::tuple<CollisionFunctor::CollisionCollectionT, bool> Simulation::collisionDetection(AutoPas_t &autopas,
                                                                                        double deltaT,
                                                                                        double collisionDistanceFactor,
                                                                                        double minDetectionRadius) {
  // pairwise interaction
  CollisionFunctor collisionFunctor(autopas.getCutoff(), deltaT, collisionDistanceFactor, minDetectionRadius);
  bool stillTuning = autopas.iteratePairwise(&collisionFunctor);
  return {collisionFunctor.getCollisions(), stillTuning};
}

size_t Simulation::simulationLoop(AutoPas_t &autopas,
                                  YoshidaIntegrator<AutoPas_t> &integrator,
                                  std::vector<Constellation> &constellations,
                                  ConfigReader &config,
                                  DomainDecomposition &domainDecomposition) {
  const auto tuningMode = config.get<bool>("autopas/tuningMode", false);
  const auto constellationInsertionFrequency = config.get<int>("io/constellationFrequency", 1);
  const auto constellationCutoff = config.get<double>("io/constellationCutoff", 0.1);
  const auto progressOutputFrequency = config.get<int>("io/progressOutputFrequency", 50);
  const auto deltaT = config.get<double>("sim/deltaT");
  const auto collisionDistanceFactor = config.get<double>("sim/collisionDistanceFactor", 1.);
  const auto timestepsPerCollisionDetection = config.get<size_t>("sim/timestepsPerCollisionDetection", 1);
  const auto decompositionType = config.get<std::string>("sim/decompositionType", "Altitude");
  // if we start from a checkpoint we want to start at the checkpoint iteration +1
  // otherwise we start at iteration 0.
  const auto startingIteration = config.getFirstIterationNr();
  const auto lastIteration = config.getLastIterationNr();
  if (timestepsPerCollisionDetection < 1) {
    SPDLOG_LOGGER_CRITICAL(
        logger.get(), "sim/timestepsPerCollisionDetection is {} but must not be <1!", timestepsPerCollisionDetection);
  }
  const auto minDetectionRadius = config.get<double>("sim/minDetectionRadius", 0.05);
  const auto minAltitude = config.get<double>("sim/minAltitude", 150.);
  std::vector<Particle> delayedInsertion;

  const auto vtkWriteFrequency = config.get<size_t>("io/vtk/writeFrequency", 0ul);

  const auto [hdf5WriteFrequency, hdf5Writer, conjuctionWriter] = initWriter(config);

  // set constellation particle IDs and fetch maxExistingParticleId
  setConstellationIDs(autopas, constellations);
  // only add the breakup model if enabled via yaml
  const std::unique_ptr<BreakupWrapper> breakupWrapper =
      config.get<bool>("sim/breakup/enabled")
          ? std::make_unique<BreakupWrapper>(config, autopas, autopas.getNumberOfParticles())
          : nullptr;

  const auto timeout = computeTimeout(config);

  size_t totalConjunctions{0ul};

  std::unique_ptr<DecompositionLogger> decompositionLogger{};

  if (decompositionType == "Altitude") {
    const auto *gridDecomp = dynamic_cast<AltitudeBasedDecomposition *>(&domainDecomposition);
  } else if (decompositionType == "RegularGrid") {
    const auto *regularGridDecomp = dynamic_cast<const RegularGridDecomposition *>(&domainDecomposition);
    decompositionLogger = std::make_unique<RegularGridDecompositionLogger>(config, *regularGridDecomp);
  } else {
    throw std::runtime_error("Unknown decomposition type: " + decompositionType);
  }

  // status output at start of the simulation
  config.printParsedValues();
  printNumParticlesPerRank(autopas, domainDecomposition);

  // in tuning mode ignore the iteration counter
  for (size_t iteration = startingIteration; iteration <= lastIteration or tuningMode; ++iteration) {
    // in case we have completed some iterations without a collisions
    // we set parentIDs to max to disable spawn protection for particles
    // recently created through breakups
    if (iterationsSinceLastCollision == 100) {
      removeParticleSpawnProtection(autopas);
    }

    // update positions
    timers.integrator.start();
    integrator.integrate(false);
    timers.integrator.stop();

    // everything below a threshold now burns up
    timers.burnUps.start();
    deleteBurnUps(autopas, minAltitude);
    timers.burnUps.stop();

    timers.constellationInsertion.start();
    // new satellites from constellations inserted over time
    if (iteration % constellationInsertionFrequency == 0) {
      updateConstellation(autopas, constellations, delayedInsertion, constellationCutoff);
    }
    timers.constellationInsertion.stop();

    timers.containerUpdate.start();
    // (potentially) update the internal data structure and collect particles which are leaving the container.

    auto leavingParticles = autopas.updateContainer();
    // for altitude decomp we cannot rely on autopas square boxes, thus
    // we need to update the leaving particles manually
    if (decompositionType == "Altitude") {
      leavingParticles = domainDecomposition.getLeavingParticles(autopas);
    }
    timers.containerUpdate.stop();

    timers.particleCommunication.start();
    auto incomingParticles = RankMigration::communicateParticles(leavingParticles, autopas, domainDecomposition);
    timers.particleCommunication.stop();

    timers.collisionDetectionImmigrants.start();
    auto collisions = RankMigration::collisionDetectionImmigrants(autopas,
                                                                  incomingParticles,
                                                                  deltaT * autopas.getVerletRebuildFrequency(),
                                                                  8.,
                                                                  collisionDistanceFactor,
                                                                  minDetectionRadius);

    timers.collisionDetectionImmigrants.stop();
    totalConjunctions += collisions.size();

    if (breakupWrapper) {
      // all particles which are part of a collision will be deleted in the breakup.
      // As we pass particle-pointers to the vector and not to autopas we have to mark the particles as deleted
      // manually
      for (const auto &[p1, p2, _, __] : collisions) {
        p1->setOwnershipState(autopas::OwnershipState::dummy);
        p2->setOwnershipState(autopas::OwnershipState::dummy);
      }
    }
    // all particles still need to be added, even if marked as deleted to not mess up autopas' internal counters.
    for (const auto &p : incomingParticles) {
      // std::cout << "Adding particle {}" << p.toString() << std::endl;
      autopas.addParticle(p);
    }

    processCollisions(iteration, collisions, *conjuctionWriter, breakupWrapper.get());

    // sanity check: after communication there should be no leaving particles left
    if (not leavingParticles.empty()) {
      SPDLOG_LOGGER_ERROR(logger.get(),
                          "It {} | {} Particles are escaping! \n{}",
                          iteration,
                          leavingParticles.size(),
                          leavingParticles);
    }

    if (iteration % timestepsPerCollisionDetection == 0) {
      timers.collisionDetection.start();
      auto [collisions, stillTuning] = collisionDetection(autopas,
                                                          deltaT * static_cast<double>(timestepsPerCollisionDetection),
                                                          collisionDistanceFactor,
                                                          minDetectionRadius);
      timers.collisionDetection.stop();
      totalConjunctions += collisions.size();
      if (collisions.size() > 0) {
        iterationsSinceLastCollision = 0;
      }

      if (tuningMode and not stillTuning) {
        dumpCalibratedConfig(config, autopas);
        return totalConjunctions;
      }

      processCollisions(iteration, collisions, *conjuctionWriter, breakupWrapper.get());
    }

    // check if we hit the timeout and abort the loop if necessary
    if (timeout != 0) {
      // quickly interrupt timers.total to update its internal total time.
      timers.total.stop();
      timers.total.start();
      // convert timer value from ns to s
      const size_t secondsSinceStart = timers.total.getTotalTime() / static_cast<size_t>(1e9);
      if (secondsSinceStart > timeout) {
        // FIXME: MPI SYNC THIS?
        SPDLOG_LOGGER_INFO(logger.get(),
                           "Simulation timeout hit! Time since simulation start ({} s) > Timeout ({} s))",
                           secondsSinceStart,
                           timeout);
        // set the config to the number of completed iterations (hence no +/-1) for the timer calculations.
        config.setValue("sim/iterations", iteration);
        // abort the loop by increasing the loop counter. This also leads to triggering the visualization
        iteration = lastIteration;
      }
    }

    timers.output.start();
    ++iterationsSinceLastCollision;
    if (iteration % progressOutputFrequency == 0 or iteration == lastIteration) {
      printProgressOutput(
          iteration, autopas.getNumberOfParticles(), totalConjunctions, domainDecomposition.getCommunicator());
    }
    // Visualization:
    if (vtkWriteFrequency and (iteration % vtkWriteFrequency == 0 or iteration == lastIteration)) {
      VTUWriter vtuWriter(config, iteration, domainDecomposition);
      // one rank has to produce the meta files
      int rank{};
      autopas::AutoPas_MPI_Comm_rank(domainDecomposition.getCommunicator(), &rank);
      if (rank == 0) {
        vtuWriter.writePVTU(config, iteration, domainDecomposition);
        decompositionLogger->writeMetafile(iteration);
      }

      vtuWriter.writeVTU(autopas);
      decompositionLogger->writePayload(iteration, autopas.getCurrentConfig());
    }
    if (hdf5WriteFrequency and (iteration % hdf5WriteFrequency == 0 or iteration == lastIteration)) {
      hdf5Writer->writeParticles(iteration, autopas);
    }
    timers.output.stop();
    // Somehow below line is always printed?
    // if (config.get<std::string>("sim/logLevel").compare("trace")) {
    //   printNumParticlesPerRank(autopas, domainDecomposition);
    // }
  }
  printNumParticlesPerRank(autopas, domainDecomposition);
  SPDLOG_LOGGER_INFO(logger.get(), "Total conjunctions: {}", totalConjunctions);
  return totalConjunctions;
}

void Simulation::removeParticleSpawnProtection(autopas::AutoPas<Particle> &autopas) {
  SPDLOG_LOGGER_INFO(logger.get(), "Removing spawn protection for all particles.");
  autopas.forEach([](auto &particle) { particle.setParentIdentifier(std::numeric_limits<size_t>::max()); });
}

std::vector<Particle> Simulation::checkedInsert(autopas::AutoPas<Particle> &autopas,
                                                const std::vector<Particle> &newSatellites,
                                                double constellationCutoff) {
  std::vector<Particle> delayedInsertion = {};

  const double collisionRadiusSquared = constellationCutoff * constellationCutoff;
  const std::array<double, 3> boxSpan = {constellationCutoff, constellationCutoff, constellationCutoff};
  for (const auto &satellite : newSatellites) {
    // only insert satellites, if they have a reasonable distance to other satellites
    bool collisionFree = true;
    // check a box around the insertion location
    const auto lowCorner = autopas::utils::ArrayMath::sub(satellite.getPosition(), boxSpan);
    const auto highCorner = autopas::utils::ArrayMath::add(satellite.getPosition(), boxSpan);
    for (auto iter = autopas.getRegionIterator(lowCorner, highCorner); iter.isValid() and collisionFree; ++iter) {
      const auto diff = autopas::utils::ArrayMath::sub(satellite.getPosition(), iter->getPosition());
      collisionFree = autopas::utils::ArrayMath::dot(diff, diff) > collisionRadiusSquared;
      if (not collisionFree) {
        SPDLOG_LOGGER_DEBUG(logger.get(),
                            "Satellite insertion delayed because it was too close to another!\n"
                            "Satellite to insert: {}\n"
                            "Colliding satellite: {}",
                            satellite,
                            *iter);
        delayedInsertion.push_back(satellite);
        break;
      }
    }
    if (collisionFree) {
      autopas.addParticle(satellite);
    }
  }
  return delayedInsertion;
}

void Simulation::setConstellationIDs(autopas::AutoPas<Particle> &autopas, std::vector<Constellation> &constellations) {
  int numRanks{};
  // AUTOPAS_MPI_COMM_WORLD should have the same size as the one stored in the decomposition
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numRanks);

  const auto lengthIDRange = std::numeric_limits<HDF5Definitions::IntType>::max() / (numRanks + constellations.size());

  // distribute globally unique ids for constellation satellites
  for (size_t i = 0; i < constellations.size(); ++i) {
    constellations[i].moveConstellationIDs(lengthIDRange * (numRanks + i));
  }
}

void Simulation::run(ConfigReader &config) {
  timers.total.start();

  timers.initialization.start();
  const auto decompositionType = config.get<std::string>("sim/decompositionType", "Altitude");
  std::unique_ptr<DomainDecomposition> domainDecomp{};
  if (decompositionType == "Altitude") {
    domainDecomp = std::make_unique<AltitudeBasedDecomposition>(config);
  } else if (decompositionType == "RegularGrid") {
    domainDecomp = std::make_unique<RegularGridDecomposition>(config);
  } else {
    throw std::runtime_error("Unknown decomposition type: " + decompositionType);
  }

  auto autopas = initAutoPas(config, *domainDecomp);
  // need to keep csvWriter and accumulator alive bc integrator relies on pointers to them but does not take ownership
  auto [csvWriter, accumulator, integrator] = initIntegrator(*autopas, config);
  SatelliteLoader::loadSatellites(*autopas, config, *domainDecomp);
  auto constellations = SatelliteLoader::loadConstellations(config, logger);
  timers.initialization.stop();

  timers.simulation.start();
  simulationLoop(*autopas, *integrator, constellations, config, *domainDecomp);
  timers.simulation.stop();

  timers.total.stop();
  timers.printTimers(config, *domainDecomp);
}

void Simulation::dumpCalibratedConfig(ConfigReader &config, const AutoPas_t &autopas) const {
  auto autopasConfig = autopas.getCurrentConfig();
  config.setValue("autopas/Newton3", autopasConfig.newton3.to_string());
  config.setValue("autopas/DataLayout", autopasConfig.dataLayout.to_string());
  config.setValue("autopas/Container", autopasConfig.container.to_string());
  config.setValue("autopas/Traversal", autopasConfig.traversal.to_string());
  config.setValue("autopas/tuningMode", "off");
  config.dumpConfig("calibratedConfig.yaml");
}

void Simulation::deleteBurnUps(autopas::AutoPas<Particle> &autopas, double burnUpAltitude) const {
  // helper function for applying comparison operators to arrays (taking OR of results)
  // => True if op is true for at least one comparison
  auto arrayComp = [](const auto &a, const auto &b, auto op) {
    bool result = false;
    for (size_t i = 0; i < a.size(); ++i) {
      result |= op(a[i], b[i]);
    }
    return result;
  };

  // skip if local box does not contain burn-up zone
  const auto critAltitude = burnUpAltitude + Physics::R_EARTH;
  const std::array<double, 3> critAltitudeBoxMin = {-critAltitude, -critAltitude, -critAltitude};
  const std::array<double, 3> critAltitudeBoxMax = {critAltitude, critAltitude, critAltitude};
  if (arrayComp(autopas.getBoxMax(), critAltitudeBoxMin, std::less()) or
      arrayComp(autopas.getBoxMin(), critAltitudeBoxMax, std::greater())) {
    return;
  }

  const auto critAltitudeSquared = critAltitude * critAltitude;
#pragma omp parallel
  for (auto particleIter = autopas.getRegionIterator({-critAltitude, -critAltitude, -critAltitude},
                                                     {critAltitude, critAltitude, critAltitude});
       particleIter != autopas.end();
       ++particleIter) {
    // Altitude above earth core
    const auto &pos = particleIter->getPosition();
    const auto particleAltitudeSquared = autopas::utils::ArrayMath::dot(pos, pos);
    if (particleAltitudeSquared < critAltitudeSquared) {
      autopas.deleteParticle(particleIter);
      SPDLOG_LOGGER_DEBUG(
          logger.get(), "Particle too close to the ground. Considered to be burning up!\n{}", *particleIter);
    }
  }
}

size_t Simulation::computeTimeout(ConfigReader &config) {
  // parse and directly convert to seconds
  constexpr bool suppressWarnings = true;
  const auto seconds = config.get<size_t>("sim/timeout/seconds", 0, suppressWarnings);
  const auto minutes = static_cast<size_t>(config.get<double>("sim/timeout/minutes", 0., suppressWarnings) * 60.);
  const auto hours = static_cast<size_t>(config.get<double>("sim/timeout/hours", 0., suppressWarnings) * 60. * 60.);
  const auto days = static_cast<size_t>(config.get<double>("sim/timeout/days", 0., suppressWarnings) * 24. * 60 * 60.);
  const auto sum = seconds + minutes + hours + days;
  // if everything resolved to 0 return the error value
  return sum;
}

void Simulation::processCollisions(size_t iteration,
                                   const CollisionFunctor::CollisionCollectionT &collisions,
                                   ConjuctionWriterInterface &conjunctionWriter,
                                   BreakupWrapper *breakupWrapper) {
  std::cout << "Processing collisions in " << iteration << std::endl;
  if (not collisions.empty()) {
    std::cout << collisions.size() << " collisions found" << std::endl;
    SPDLOG_LOGGER_TRACE(logger.get(), "The following particles collided between ranks:");
    for (const auto &[p1, p2, _, __] : collisions) {
      std::cout << "Collision between " << p1->toString() << " and " << p2->toString() << std::endl;
      SPDLOG_LOGGER_TRACE(logger.get(), "({}, {})", p1->toString(), p2->toString());
    }
    iterationsSinceLastCollision = 0;
  }
  timers.collisionWriting.start();
  conjunctionWriter.writeConjunctions(iteration, collisions);
  timers.collisionWriting.stop();

  if (breakupWrapper and not collisions.empty()) {
    timers.collisionSimulation.start();
    breakupWrapper->simulateBreakup(collisions);
    timers.collisionSimulation.stop();
  }
}

void Simulation::printNumParticlesPerRank(const autopas::AutoPas<Particle> &autopas,
                                          const DomainDecomposition &decomposition) const {
  const auto communicator = decomposition.getCommunicator();
  int myRank{};
  autopas::AutoPas_MPI_Comm_rank(communicator, &myRank);
  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(communicator, &numRanks);
  std::vector<size_t> numParticlesPerRank(numRanks);
  // rank zero collects all data and prints
  if (myRank == 0) {
    numParticlesPerRank[myRank] = autopas.getNumberOfParticles();
    for (int rank = 1; rank < numRanks; ++rank) {
      autopas::AutoPas_MPI_Recv(&numParticlesPerRank[rank],
                                1,
                                AUTOPAS_MPI_UNSIGNED_LONG,
                                rank,
                                rank,
                                communicator,
                                AUTOPAS_MPI_STATUS_IGNORE);
    }
    const auto numParticlesAvg =
        std::accumulate(numParticlesPerRank.begin(), numParticlesPerRank.end(), 0ul) / numRanks;
    logger.log(LADDS::Logger::Level::info,
               "Particles per rank (Avg: {}) : {}",
               numParticlesAvg,
               autopas::utils::ArrayUtils::to_string(numParticlesPerRank, " ", {"", ""}));
  } else {
    // other ranks only send
    const auto myNumParticles = autopas.getNumberOfParticles();
    autopas::AutoPas_MPI_Send(&myNumParticles, 1, AUTOPAS_MPI_UNSIGNED_LONG, 0, myRank, communicator);
  }
}

void Simulation::printProgressOutput(size_t iteration,
                                     size_t numParticlesLocal,
                                     size_t totalConjunctionsLocal,
                                     autopas::AutoPas_MPI_Comm const &comm) {
  unsigned long numParticlesGlobal{};
  unsigned long totalConjunctionsGlobal{};
  autopas::AutoPas_MPI_Reduce(
      &numParticlesLocal, &numParticlesGlobal, 1, AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_SUM, 0, comm);
  autopas::AutoPas_MPI_Reduce(
      &totalConjunctionsLocal, &totalConjunctionsGlobal, 1, AUTOPAS_MPI_UNSIGNED_LONG, AUTOPAS_MPI_SUM, 0, comm);
  SPDLOG_LOGGER_INFO(logger.get(),
                     "It {} | Total particles: {} | Total conjunctions: {}",
                     iteration,
                     numParticlesGlobal,
                     totalConjunctionsGlobal);
}

}  // namespace LADDS