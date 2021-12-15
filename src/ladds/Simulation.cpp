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
void setAutoPasOption(const YAML::Node &node, F setterFun, const std::set<Option> &defaultVal) {
  if (node.IsDefined()) {
    const auto options = Option::parseOptions(node.as<std::string>());
    setterFun(options);
  } else {
    setterFun(defaultVal);
  }
}
}  // namespace

std::unique_ptr<AutoPas_t> Simulation::initAutoPas(const YAML::Node &config) {
  auto autopas = std::make_unique<AutoPas_t>();

  const auto maxAltitude = config["sim"]["maxAltitude"].as<double>();
  const auto cutoff = config["autopas"]["cutoff"].as<double>();
  const auto verletSkin = config["autopas"]["skin"].as<double>();
  const auto verletRebuildFrequency = config["autopas"]["rebuildFrequency"].as<unsigned int>();
  const auto desiredCellsPerDimension = config["autopas"]["desiredCellsPerDimension"].as<double>();

  autopas->setBoxMin({-maxAltitude, -maxAltitude, -maxAltitude});
  autopas->setBoxMax({maxAltitude, maxAltitude, maxAltitude});
  autopas->setCutoff(cutoff);
  autopas->setVerletSkin(verletSkin);
  autopas->setVerletRebuildFrequency(verletRebuildFrequency);
  // Scale Cell size so that we get the desired number of cells
  // -2 because internally there will be two halo cells added on top of maxAltitude
  autopas->setCellSizeFactor((maxAltitude * 2.) / ((cutoff + verletSkin) * (desiredCellsPerDimension - 2)));

  // only restrict AutoPas options if we are not in tuning mode
  const auto tuningModeNode = config["autopas"]["tuningMode"];
  const auto tuningMode = tuningModeNode.IsDefined() and tuningModeNode.as<bool>();
  if (not tuningMode) {
    setAutoPasOption<autopas::Newton3Option>(config["autopas"]["Newton3"],
                                             [&](const auto &op) { autopas->setAllowedNewton3Options(op); },
                                             {autopas::Newton3Option::enabled});
    setAutoPasOption<autopas::DataLayoutOption>(config["autopas"]["DataLayout"],
                                                [&](const auto &op) { autopas->setAllowedDataLayouts(op); },
                                                {autopas::DataLayoutOption::aos});
    setAutoPasOption<autopas::ContainerOption>(config["autopas"]["Container"],
                                               [&](const auto &op) { autopas->setAllowedContainers(op); },
                                               {autopas::ContainerOption::varVerletListsAsBuild});
    setAutoPasOption<autopas::TraversalOption>(config["autopas"]["Traversal"],
                                               [&](const auto &op) { autopas->setAllowedTraversals(op); },
                                               {autopas::TraversalOption::vvl_as_built});
  }
  // arbitrary number. Can be changed to whatever makes sense.
  autopas->setTuningInterval(10000);
  autopas->setSelectorStrategy(autopas::SelectorStrategyOption::fastestMean);
  autopas->setNumSamples(verletRebuildFrequency);
  if (config["autopas"]["logLevel"].IsDefined()) {
    autopas::Logger::get()->set_level(spdlog::level::from_str(config["autopas"]["logLevel"].as<std::string>()));
  }
  autopas->init();

  return autopas;
}

std::tuple<std::unique_ptr<FileOutput<AutoPas_t>>,
           std::unique_ptr<Acceleration::AccelerationAccumulator<AutoPas_t>>,
           std::unique_ptr<Integrator<AutoPas_t>>>
Simulation::initIntegrator(AutoPas_t &autopas, const YAML::Node &config) {
  // initialization of the integrator
  std::array<bool, 8> selectedPropagatorComponents{};
  const auto &configProp = config["sim"]["prop"];
  if (configProp.IsDefined()) {
    selectedPropagatorComponents = {
        configProp["useKEPComponent"].IsDefined() and configProp["useKEPComponent"].as<bool>(),
        configProp["useJ2Component"].IsDefined() and configProp["useJ2Component"].as<bool>(),
        configProp["useC22Component"].IsDefined() and configProp["useC22Component"].as<bool>(),
        configProp["useS22Component"].IsDefined() and configProp["useS22Component"].as<bool>(),
        configProp["useSOLComponent"].IsDefined() and configProp["useSOLComponent"].as<bool>(),
        configProp["useLUNComponent"].IsDefined() and configProp["useLUNComponent"].as<bool>(),
        configProp["useSRPComponent"].IsDefined() and configProp["useSRPComponent"].as<bool>(),
        configProp["useDRAGComponent"].IsDefined() and configProp["useDRAGComponent"].as<bool>()};
  }

  auto csvWriter = std::make_unique<FileOutput<AutoPas_t>>(
      autopas, config["io"]["output_file"].as<std::string>(), OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_unique<Acceleration::AccelerationAccumulator<AutoPas_t>>(
      selectedPropagatorComponents, autopas, 0.0, *csvWriter);
  auto deltaT = config["sim"]["deltaT"].as<double>();
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

void Simulation::collisionDetection(size_t iteration,
                                    AutoPas_t &autopas,
                                    ConjunctionLogger &conjunctionLogger,
                                    size_t &totalConjunctions,
                                    size_t progressOutputFrequency,
                                    double deltaT,
                                    double conjunctionThreshold) {
  // pairwise interaction
  CollisionFunctor collisionFunctor(autopas.getCutoff(), deltaT, conjunctionThreshold);
  autopas.iteratePairwise(&collisionFunctor);
  auto collisions = collisionFunctor.getCollisions();
  for (const auto &[p1, p2AndDistanceSquare] : collisions) {
    totalConjunctions++;
    const auto &[p2, distanceSquare] = p2AndDistanceSquare;
    conjunctionLogger.log(iteration, *p1, *p2, distanceSquare);
    SPDLOG_LOGGER_DEBUG(logger.get(), "{} | {} | distanceSquare={}", p1->getID(), p2->getID(), distanceSquare);
  }
  if (iteration % progressOutputFrequency == 0) {
    SPDLOG_LOGGER_INFO(
        logger.get(), "It {} - Encounters:{} Total conjunctions:{}", iteration, collisions.size(), totalConjunctions);
  }
}

void Simulation::simulationLoop(AutoPas_t &autopas,
                                Integrator<AutoPas_t> &integrator,
                                std::vector<Constellation> &constellations,
                                const YAML::Node &config) {
  const auto iterations = config["sim"]["iterations"].as<size_t>();
  const auto vtkWriteFrequency = config["io"]["vtkWriteFrequency"].as<size_t>();
  const auto constellationInsertionFrequency =
      config["io"]["constellationFrequency"].IsDefined() ? config["io"]["constellationFrequency"].as<int>() : 1;
  const auto constellationCutoff =
      config["io"]["constellationCutoff"].IsDefined() ? config["io"]["constellationCutoff"].as<double>() : 0.1;
  const auto progressOutputFrequency =
      config["io"]["progressOutputFrequency"].IsDefined() ? config["io"]["progressOutputFrequency"].as<int>() : 50;
  const auto deltaT = config["sim"]["deltaT"].as<double>();
  const auto conjunctionThreshold = config["sim"]["conjunctionThreshold"].as<double>();
  std::vector<Particle> delayedInsertion;

  size_t totalConjunctions{0ul};
  ConjunctionLogger conjunctionLogger("");

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
    collisionDetection(
        i, autopas, conjunctionLogger, totalConjunctions, progressOutputFrequency, deltaT, conjunctionThreshold);
    timers.collisionDetection.stop();

    // TODO insert breakup model here

    timers.output.start();
    // Visualization:
    if (i % vtkWriteFrequency == 0) {
      VTUWriter::writeVTK(i, autopas);
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

void Simulation::run(const YAML::Node &config) {
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
