/**
 * @file Simulation.cpp
 * @author F. Gratl
 * @date 30.11.21
 */

#include "Simulation.h"

#include <breakupModel/output/VTKWriter.h>
#include <ladds/io/DatasetReader.h>
#include <ladds/particle/SatelliteToParticleConverter.h>

#include <array>
#include <tuple>
#include <vector>

#include "ladds/io/ConjunctionLogger.h"
#include "ladds/io/hdf5/HDF5Writer.h"
#include "ladds/io/SatelliteLoader.h"
#include "ladds/io/VTUWriter.h"
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
 * @param defaultVal
 */
template <class Option, class F>
void setAutoPasOption(ConfigReader &config, const std::string &valuePath, F setterFun, const Option &defaultVal) {
  auto valStr = config.template get<std::string>(valuePath, defaultVal.to_string(), true);
  const auto options = Option::parseOptions(valStr);
  setterFun(options);
}
}  // namespace

std::unique_ptr<AutoPas_t> Simulation::initAutoPas(ConfigReader &config) {
  auto autopas = std::make_unique<AutoPas_t>();

  const auto maxAltitude = config.get<double>("sim/maxAltitude");
  const auto cutoff = config.get<double>("autopas/cutoff");
  const auto desiredCellsPerDimension = config.get<double>("autopas/desiredCellsPerDimension");
  const auto deltaT = config.get<double>("sim/deltaT");
  const auto verletRebuildFrequency = config.get<unsigned int>("autopas/rebuildFrequency");
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

  // only restrict AutoPas options if we are not in tuning mode
  const auto tuningMode = config.get<bool>("autopas/tuningMode", false);
  if (not tuningMode) {
    setAutoPasOption<autopas::Newton3Option>(config,
                                             "autopas/Newton3",
                                             [&](const auto &op) { autopas->setAllowedNewton3Options(op); },
                                             {autopas::Newton3Option::enabled});
    setAutoPasOption<autopas::DataLayoutOption>(config,
                                                "autopas/DataLayout",
                                                [&](const auto &op) { autopas->setAllowedDataLayouts(op); },
                                                {autopas::DataLayoutOption::aos});
    setAutoPasOption<autopas::ContainerOption>(config,
                                               "autopas/Container",
                                               [&](const auto &op) { autopas->setAllowedContainers(op); },
                                               {autopas::ContainerOption::varVerletListsAsBuild});
    setAutoPasOption<autopas::TraversalOption>(config,
                                               "autopas/Traversal",
                                               [&](const auto &op) { autopas->setAllowedTraversals(op); },
                                               {autopas::TraversalOption::vvl_as_built});
  }
  // arbitrary number. Can be changed to whatever makes sense.
  autopas->setTuningInterval(10000);
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
      selectedPropagatorComponents, autopas, 0.0, *csvWriter);
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

std::unordered_map<Particle *, std::tuple<Particle *, double>> Simulation::collisionDetection(
    AutoPas_t &autopas, double deltaT, double conjunctionThreshold) {
  // pairwise interaction
  CollisionFunctor collisionFunctor(autopas.getCutoff(), deltaT, conjunctionThreshold);
  autopas.iteratePairwise(&collisionFunctor);
  return collisionFunctor.getCollisions();
}

void Simulation::simulationLoop(AutoPas_t &autopas,
                                Integrator<AutoPas_t> &integrator,
                                std::vector<Constellation> &constellations,
                                ConfigReader &config) {
  const auto iterations = config.get<size_t>("sim/iterations");
  const auto constellationInsertionFrequency = config.get<int>("io/constellationFrequency", 1);
  const auto constellationCutoff = config.get<double>("io/constellationCutoff", 0.1);
  const auto progressOutputFrequency = config.get<int>("io/progressOutputFrequency", 50);
  const auto deltaT = config.get<double>("sim/deltaT");
  const auto conjunctionThreshold = config.get<double>("sim/conjunctionThreshold");
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

  for (size_t i = 0ul; i < iterations; ++i) {
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
    auto collisions = collisionDetection(autopas, deltaT, conjunctionThreshold);
    timers.collisionDetection.stop();

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
