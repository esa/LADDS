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
#include <tuple>
#include <vector>

#include "ladds/io/ConjunctionLogger.h"
#include "ladds/io/SatelliteLoader.h"
#include "ladds/io/VTUWriter.h"
#include "ladds/io/hdf5/HDF5Writer.h"
#include "ladds/particle/Constellation.h"

// Declare the main AutoPas class as extern template instantiation. It is instantiated in AutoPasClass.cpp.
extern template class autopas::AutoPas<Particle>;
extern template bool autopas::AutoPas<Particle>::iteratePairwise(CollisionFunctor *);

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

std::unique_ptr<AutoPas_t> Simulation::initAutoPas(ConfigReader &config) {
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

  autopas->setBoxMin({-maxAltitude, -maxAltitude, -maxAltitude});
  autopas->setBoxMax({maxAltitude, maxAltitude, maxAltitude});
  autopas->setCutoff(cutoff);
  autopas->setVerletSkin(verletSkin);
  autopas->setVerletRebuildFrequency(verletRebuildFrequency);
  // Scale Cell size so that we get the desired number of cells
  // -2 because internally there will be two halo cells added on top of maxAltitude
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
  autopas::Logger::get()->set_level(spdlog::level::from_str(config.get<std::string>("autopas/logLevel", "off")));

  autopas->init();

  return autopas;
}

std::tuple<std::unique_ptr<FileOutput<AutoPas_t>>,
           std::unique_ptr<Acceleration::AccelerationAccumulator<AutoPas_t>>,
           std::unique_ptr<Integrator<AutoPas_t>>>
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

  auto csvWriter = std::make_unique<FileOutput<AutoPas_t>>(autopas,
                                                           config.get<std::string>("io/output_file", "propagator.csv"),
                                                           OutputFile::CSV,
                                                           selectedPropagatorComponents);
  auto accumulator = std::make_unique<Acceleration::AccelerationAccumulator<AutoPas_t>>(
      selectedPropagatorComponents, autopas, 0.0, csvWriter.get());
  auto deltaT = config.get<double>("sim/deltaT");
  auto integrator = std::make_unique<Integrator<AutoPas_t>>(autopas, *accumulator, deltaT);

  return std::make_tuple<>(std::move(csvWriter), std::move(accumulator), std::move(integrator));
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
                                                                                        double conjunctionThreshold,
                                                                                        double minDetectionRadius) {
  // pairwise interaction
  CollisionFunctor collisionFunctor(autopas.getCutoff(), deltaT, conjunctionThreshold, minDetectionRadius);
  bool stillTuning = autopas.iteratePairwise(&collisionFunctor);
  return {collisionFunctor.getCollisions(), stillTuning};
}

void Simulation::simulationLoop(AutoPas_t &autopas,
                                Integrator<AutoPas_t> &integrator,
                                std::vector<Constellation> &constellations,
                                ConfigReader &config) {
  const auto tuningMode = config.get<bool>("autopas/tuningMode", false);
  const auto iterations = config.get<size_t>("sim/iterations");
  const auto constellationInsertionFrequency = config.get<int>("io/constellationFrequency", 1);
  const auto constellationCutoff = config.get<double>("io/constellationCutoff", 0.1);
  const auto progressOutputFrequency = config.get<int>("io/progressOutputFrequency", 50);
  const auto deltaT = config.get<double>("sim/deltaT");
  const auto conjunctionThreshold = config.get<double>("sim/conjunctionThreshold");
  const auto minDetectionRadius = config.get<double>("sim/minDetectionRadius", 0.1);
  std::vector<Particle> delayedInsertion;

  const auto vtkWriteFrequency = config.get<size_t>("io/vtk/writeFrequency", 0ul);

  auto hdf5WriteFrequency = 0u;

  std::shared_ptr<HDF5Writer> hdf5Writer;
  std::shared_ptr<ConjuctionWriterInterface> conjuctionWriter;
  if (config.defines("io/hdf5")) {
    hdf5WriteFrequency = config.get<unsigned int>("io/hdf5/writeFrequency", 0u);
    const auto hdf5CompressionLvl = config.get<unsigned int>("io/hdf5/compressionLevel", 0u);
    const auto hdf5FileName = config.get<std::string>("io/hdf5/fileName", "simulationData.h5");

    hdf5Writer = std::make_shared<HDF5Writer>(hdf5FileName, hdf5CompressionLvl);
    conjuctionWriter = hdf5Writer;
  } else {
    conjuctionWriter = std::make_shared<ConjunctionLogger>("");
  }

  size_t totalConjunctions{0ul};

  config.printParsedValues();

  // in tuning mode ignore the iteration counter
  for (size_t i = 0ul; i < iterations or tuningMode; ++i) {
    // update positions
    timers.integrator.start();
    integrator.integrate(false);
    timers.integrator.stop();

    timers.constellationInsertion.start();
    // new satellites from constellations inserted over time
    if (i % constellationInsertionFrequency == 0) {
      updateConstellation(autopas, constellations, delayedInsertion, constellationCutoff);
    }
    timers.constellationInsertion.stop();

    // TODO MPI: handle particle exchange between ranks
    timers.containerUpdate.start();
    // (potentially) update the internal data structure and collect particles which are leaving the container.
    const auto escapedParticles = autopas.updateContainer();
    timers.containerUpdate.stop();

    // sanity check
    if (not escapedParticles.empty()) {
      SPDLOG_LOGGER_ERROR(logger.get(), "Particles are escaping! \n{}", escapedParticles);
    }
    // TODO Check for particles that burn up

    timers.collisionDetection.start();
    auto [collisions, stillTuning] = collisionDetection(autopas, deltaT, conjunctionThreshold, minDetectionRadius);
    timers.collisionDetection.stop();

    if (tuningMode and not stillTuning) {
      dumpCalibratedConfig(config, autopas);
      return;
    }

    timers.collisionWriting.start();
    totalConjunctions += collisions.size();
    conjuctionWriter->writeConjunctions(i, collisions);
    if (i % progressOutputFrequency == 0) {
      SPDLOG_LOGGER_INFO(
          logger.get(), "It {} - Encounters:{} Total conjunctions:{}", i, collisions.size(), totalConjunctions);
    }
    timers.collisionWriting.stop();

    // TODO insert breakup model here

    timers.output.start();
    // Visualization:
    if (vtkWriteFrequency and i % vtkWriteFrequency == 0) {
      VTUWriter::writeVTU(i, autopas);
    }
    if (hdf5WriteFrequency and i % hdf5WriteFrequency == 0) {
      hdf5Writer->writeParticles(i, autopas);
    }
    timers.output.stop();
  }
  SPDLOG_LOGGER_INFO(logger.get(), "Total conjunctions: {}", totalConjunctions);
}

std::vector<Particle> Simulation::checkedInsert(autopas::AutoPas<Particle> &autopas,
                                                const std::vector<Particle> &newSatellites,
                                                double constellationCutoff) {
  std::vector<Particle> delayedInsertion = {};

  const double collisionRadius = constellationCutoff;
  const double collisionRadiusSquared = collisionRadius * collisionRadius;
  const std::array<double, 3> boxSpan = {collisionRadius, collisionRadius, collisionRadius};
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

void Simulation::run(ConfigReader &config) {
  timers.total.start();

  timers.initialization.start();
  auto autopas = initAutoPas(config);
  // need to keep csvWriter and accumulator alive bc integrator relies on pointers to them but does not take ownership
  auto [csvWriter, accumulator, integrator] = initIntegrator(*autopas, config);
  SatelliteLoader::loadSatellites(*autopas, config, logger);
  auto constellations = SatelliteLoader::loadConstellations(config, logger);
  timers.initialization.stop();

  timers.simulation.start();
  simulationLoop(*autopas, *integrator, constellations, config);
  timers.simulation.stop();

  timers.total.stop();
  timers.printTimers(config);
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
