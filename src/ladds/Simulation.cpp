/**
 * @file Simulation.cpp
 * @author F. Gratl
 * @date 30.11.21
 */

#include "Simulation.h"

#include <breakupModel/output/VTKWriter.h>
#include <ladds/io/DatasetReader.h>
#include <ladds/particle/SatelliteToParticleConverter.h>

#include <tuple>

// Declare the main AutoPas class as extern template instantiation. It is instantiated in AutoPasClass.cpp.
extern template class autopas::AutoPas<Particle>;
extern template bool autopas::AutoPas<Particle>::iteratePairwise(CollisionFunctor *);

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
  autopas->setAllowedNewton3Options({autopas::Newton3Option::disabled});
  autopas->setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autopas->setAllowedContainers({autopas::ContainerOption::varVerletListsAsBuild});
  autopas->setAllowedTraversals({autopas::TraversalOption::vvl_as_built});
  autopas::Logger::get()->set_level(spdlog::level::from_str(config["autopas"]["logLevel"].as<std::string>()));
  autopas->init();

  return autopas;
}

auto Simulation::initIntegrator(AutoPas_t &autopas, const YAML::Node &config) {
  // initialization of the integrator
  std::array<bool, 8> selectedPropagatorComponents{config["sim"]["prop"]["useKEPComponent"].as<bool>(),
                                                   config["sim"]["prop"]["useJ2Component"].as<bool>(),
                                                   config["sim"]["prop"]["useC22Component"].as<bool>(),
                                                   config["sim"]["prop"]["useS22Component"].as<bool>(),
                                                   config["sim"]["prop"]["useSOLComponent"].as<bool>(),
                                                   config["sim"]["prop"]["useLUNComponent"].as<bool>(),
                                                   config["sim"]["prop"]["useSRPComponent"].as<bool>(),
                                                   config["sim"]["prop"]["useDRAGComponent"].as<bool>()};

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
  auto actualSatellites = DatasetReader::readDataset(config["io"]["posFileName"].as<std::string>(),
                                                     config["io"]["velFileName"].as<std::string>());
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

void Simulation::simulationLoop(AutoPas_t &autopas, Integrator<AutoPas_t> &integrator, const YAML::Node &config) {
  const auto cutoff = config["autopas"]["cutoff"].as<double>();
  const auto iterations = config["sim"]["iterations"].as<size_t>();
  const auto vtkWriteFrequency = config["io"]["vtkWriteFrequency"].as<size_t>();

  for (size_t i = 0; i < iterations; ++i) {
    // update positions
    timers.integrator.start();
    integrator.integrate(false);
    timers.integrator.stop();

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
    CollisionFunctor collisionFunctor(cutoff);
    autopas.iteratePairwise(&collisionFunctor);
    auto collisions = collisionFunctor.getCollisions();
    SPDLOG_LOGGER_INFO(logger.get(), "Iteration {} - Close encounters: {}", i, collisions.size());
    for (const auto &[p1, p2] : collisions) {
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

void Simulation::printTimers(const YAML::Node &config) const {
  const auto iterations = config["sim"]["iterations"].as<size_t>();

  const auto timeTotal = timers.total.getTotalTime();
  const auto timeSim = timers.simulation.getTotalTime();
  const auto maximumNumberOfDigits = static_cast<int>(std::to_string(timeTotal).length());
  std::cout << timerToString("Total                   ", timeTotal, maximumNumberOfDigits);
  std::cout << timerToString(
      "  Initialization        ", timers.initialization.getTotalTime(), maximumNumberOfDigits, timeTotal);
  std::cout << timerToString("  Simulation            ", timeSim, maximumNumberOfDigits, timeTotal);
  std::cout << timerToString(
      "    Integrator          ", timers.integrator.getTotalTime(), maximumNumberOfDigits, timeSim);
  std::cout << timerToString(
      "    Collision detection ", timers.collisionDetection.getTotalTime(), maximumNumberOfDigits, timeSim);
  std::cout << timerToString(
      "    Container update    ", timers.containerUpdate.getTotalTime(), maximumNumberOfDigits, timeSim);
  std::cout << timerToString(
      "    Output              ", timers.output.getTotalTime(), maximumNumberOfDigits, timeTotal);
  std::cout << timerToString("One iteration           ", timeSim / iterations, maximumNumberOfDigits, timeTotal);
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
  timers.initialization.stop();

  timers.simulation.start();
  simulationLoop(*autopas, *integrator, config);
  timers.simulation.stop();

  timers.total.stop();
  printTimers(config);
}