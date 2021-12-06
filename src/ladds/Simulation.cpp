/**
 * @file Simulation.cpp
 * @author F. Gratl
 * @date 30.11.21
 */

#include <breakupModel/output/VTKWriter.h>
#include <ladds/io/DatasetReader.h>
#include <ladds/particle/SatelliteToParticleConverter.h>

#include <array>
#include <tuple>
#include <vector>

#include "Simulation.h"
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

std::unique_ptr<Simulation::AutoPas_t> Simulation::initAutoPas(const YAML::Node &config) {
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

std::tuple<std::unique_ptr<FileOutput<Simulation::AutoPas_t>>,
           std::unique_ptr<Acceleration::AccelerationAccumulator<Simulation::AutoPas_t>>,
           std::unique_ptr<Integrator<Simulation::AutoPas_t>>>
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

void Simulation::loadSatellites(AutoPas_t &autopas, const YAML::Node &config) {
  // Read in scenario
  auto actualSatellites =
      DatasetReader::readDataset(std::string(DATADIR) + config["io"]["posFileName"].as<std::string>(),
                                 std::string(DATADIR) + config["io"]["velFileName"].as<std::string>());
  SPDLOG_LOGGER_DEBUG(logger.get(), "Parsed {} satellites", actualSatellites.size());

  const auto maxAltitude = config["sim"]["maxAltitude"].as<double>();
  double minAltitudeFound{std::numeric_limits<double>::max()};
  double maxAltitudeFound{0.};
  // Convert satellites to particles
  for (const auto &particle : actualSatellites) {
    auto pos = particle.getPosition();
    double altitude = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
    minAltitudeFound = std::min(minAltitudeFound, altitude);
    maxAltitudeFound = std::max(maxAltitudeFound, altitude);
    if (altitude < maxAltitude) {
      autopas.addParticle(particle);
    }
  }
  SPDLOG_LOGGER_INFO(logger.get(), "Min altitude is {}", minAltitudeFound);
  SPDLOG_LOGGER_INFO(logger.get(), "Max altitude is {}", maxAltitudeFound);
  SPDLOG_LOGGER_INFO(logger.get(), "Number of particles: {}", autopas.getNumberOfParticles());
}

std::vector<Constellation> Simulation::loadConstellations(const YAML::Node &config) {
  std::vector<Constellation> constellations;

  if (config["io"]["constellationList"].IsDefined()) {
    const auto insertionFrequency =
        config["io"]["constellationFrequency"].IsDefined() ? config["io"]["constellationFrequency"].as<int>() : 1;
    auto constellationDataStr = config["io"]["constellationList"].as<std::string>();
    // count constellation by counting ';'
    int nConstellations = 1;
    for (char con : constellationDataStr) {
      if (con == ';') {
        nConstellations++;
      }
    }

    // parse constellation info
    constellations.reserve(nConstellations);
    for (int i = 0; i < nConstellations; ++i) {
      unsigned long offset = constellationDataStr.find(';', 0);
      if (offset == 0) {
        constellations.emplace_back(Constellation(constellationDataStr, insertionFrequency));
        break;
      } else {
        constellations.emplace_back(Constellation(constellationDataStr.substr(0, offset), insertionFrequency));
        constellationDataStr.erase(0, offset + 1);
      }
    }

    size_t constellationTotalNumSatellites = 0;
    for (const auto &constellation : constellations) {
      constellationTotalNumSatellites += constellation.getConstellationSize();
    }

    SPDLOG_LOGGER_INFO(logger.get(),
                       "{} more particles will be added from {} constellations",
                       constellationTotalNumSatellites,
                       nConstellations);
  }
  return constellations;
}

void Simulation::simulationLoop(AutoPas_t &autopas,
                                Integrator<AutoPas_t> &integrator,
                                std::vector<Constellation> &constellations,
                                const YAML::Node &config) {
  const auto cutoff = config["autopas"]["cutoff"].as<double>();
  const auto iterations = config["sim"]["iterations"].as<size_t>();
  const auto vtkWriteFrequency = config["io"]["vtkWriteFrequency"].as<size_t>();
  const auto constellationInsertionFrequency =
      config["io"]["constellationFrequency"].IsDefined() ? config["io"]["constellationFrequency"].as<int>() : 1;
  std::vector<Particle> delayedInsertion;

  for (size_t i = 0ul; i < iterations; ++i) {
    // update positions
    timers.integrator.start();
    integrator.integrate(false);
    timers.integrator.stop();

    timers.constellationInsertion.start();
    // new satellites from constellations inserted over time
    if (i % constellationInsertionFrequency == 0) {
      for (auto &constellation : constellations) {
        // new satellites are gradually added to the simulation according to their starting time and operation duration
        auto newSatellites = constellation.tick();

        // add waiting satellites to newSatellites
        newSatellites.insert(newSatellites.end(), delayedInsertion.begin(), delayedInsertion.end());
        delayedInsertion = checkedInsert(autopas, newSatellites, cutoff);
      }
    }
    timers.constellationInsertion.stop();

    // TODO MPI: handle particle exchange between ranks
    timers.containerUpdate.start();
    // (potentially) update the internal data structure and collect particles which are leaving the container.
    const auto escapedParticles = autopas.updateContainer();
    timers.containerUpdate.stop();

    if (not escapedParticles.empty()) {
      SPDLOG_LOGGER_ERROR(logger.get(), "Particles are escaping! \n{}", escapedParticles);
    }
    // TODO Check for particles that burn up

    timers.collisionDetection.start();
    // pairwise interaction
    CollisionFunctor collisionFunctor(cutoff, config["sim"]["deltaT"], 0.1 * cutoff);
    autopas.iteratePairwise(&collisionFunctor);
    auto collisions = collisionFunctor.getCollisions();
    SPDLOG_LOGGER_INFO(logger.get(), "Iteration {} - Close encounters: {}", i, collisions.size());
    for (const auto &[p1, p2AndDist] : collisions) {
      const auto &[p2, dist] = p2AndDist;
      SPDLOG_LOGGER_DEBUG(logger.get(), "{} | {}", p1->getID(), p2->getID());
    }
    timers.collisionDetection.stop();

    // TODO insert breakup model here

    timers.output.start();
    // Visualization:
    if (i % vtkWriteFrequency == 0) {
      VTKWriter vtkWriter("output_" + std::to_string(i) + ".vtu");
      std::vector<Satellite> allParticles;
      allParticles.reserve(autopas.getNumberOfParticles());
      for (const auto &p : autopas) {
        allParticles.push_back(SatelliteToParticleConverter::convertParticleToSatellite(p));
      }
      // sort particles by Id to provide consistent output files
      std::sort(allParticles.begin(), allParticles.end(), [](const auto &p1, const auto &p2) {
        return p1.getId() < p2.getId();
      });
      vtkWriter.printResult(allParticles);
    }
    timers.output.stop();
  }
}

std::vector<Particle> Simulation::checkedInsert(autopas::AutoPas<Particle> &autopas,
                                                const std::vector<Particle> &newSatellites,
                                                double cutoff) {
  std::vector<Particle> delayedInsertion = {};

  const double collisionRadius = 2 * cutoff;
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

void Simulation::printTimers(const YAML::Node &config) const {
  const auto iterations = config["sim"]["iterations"].as<size_t>();

  const auto timeTotal = timers.total.getTotalTime();
  const auto timeSim = timers.simulation.getTotalTime();
  const auto maximumNumberOfDigits = static_cast<int>(std::to_string(timeTotal).length());
  std::cout << timerToString("Total                       ", timeTotal, maximumNumberOfDigits);
  std::cout << timerToString(
      "  Initialization            ", timers.initialization.getTotalTime(), maximumNumberOfDigits, timeTotal);
  std::cout << timerToString("  Simulation            ", timeSim, maximumNumberOfDigits, timeTotal);
  std::cout << timerToString(
      "    Integrator              ", timers.integrator.getTotalTime(), maximumNumberOfDigits, timeSim);
  std::cout << timerToString(
      "    Constellation insertion ", timers.constellationInsertion.getTotalTime(), maximumNumberOfDigits, timeSim);
  std::cout << timerToString(
      "    Collision detection     ", timers.collisionDetection.getTotalTime(), maximumNumberOfDigits, timeSim);
  std::cout << timerToString(
      "    Container update        ", timers.containerUpdate.getTotalTime(), maximumNumberOfDigits, timeSim);
  std::cout << timerToString(
      "    Output                  ", timers.output.getTotalTime(), maximumNumberOfDigits, timeTotal);
  std::cout << timerToString("One iteration               ", timeSim / iterations, maximumNumberOfDigits, timeTotal);
}

/**
 * Turns the timers into a human readable string.
 * @param name: The timer's name.
 * @param timeNS: The time in nano seconds.
 * @param numberWidth: The precision of the printed number.
 * @param maxTime: The simulation's total execution time.
 * @return All information of the timer in a human readable string.
 *
 * @note Taken from md-flexible.
 */
std::string Simulation::timerToString(const std::string &name, long timeNS, int numberWidth, long maxTime) {
  // only print timers that were actually used
  if (timeNS == 0) {
    return "";
  }

  std::ostringstream ss;
  ss << std::fixed << std::setprecision(floatStringPrecision) << name << " : " << std::setw(numberWidth) << std::right
     << timeNS
     << " ns ("
     // min width of the representation of seconds is numberWidth - 9 (from conversion) + 4 (for dot and digits after)
     << std::setw(numberWidth - 5) << ((double)timeNS * 1e-9) << "s)";
  if (maxTime != 0) {
    ss << " =" << std::setw(7) << std::right << ((double)timeNS / (double)maxTime * 100) << "%";
  }
  ss << std::endl;
  return ss.str();
}

void Simulation::run(const YAML::Node &config) {
  timers.total.start();

  timers.initialization.start();
  auto autopas = initAutoPas(config);
  // need to keep csvWriter and accumulator alive bc integrator relies on pointers to them but does not take ownership
  auto [csvWriter, accumulator, integrator] = initIntegrator(*autopas, config);
  loadSatellites(*autopas, config);
  auto constellations = loadConstellations(config);
  timers.initialization.stop();

  timers.simulation.start();
  simulationLoop(*autopas, *integrator, constellations, config);
  timers.simulation.stop();

  timers.total.stop();
  printTimers(config);
}