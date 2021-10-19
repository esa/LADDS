/**
 * @file CollisionFunctor.h
 * @author F. Gratl
 * @date 28.06.21
 */

#include <autopas/AutoPasDecl.h>
#include <autopas/utils/Timer.h>
#include <breakupModel/output/VTKWriter.h>
#include <satellitePropagator/io/FileOutput.h>
#include <satellitePropagator/physics/AccelerationAccumulator.h>
#include <satellitePropagator/physics/Integrator.h>

#include <iostream>

#include "CollisionFunctor.h"
#include "DatasetReader.h"
#include "LoadConfig.h"
#include "Logger.h"
#include "Particle.h"
#include "SatelliteToParticleConverter.h"
#include "spdlog/fmt/ostr.h"

// Declare the main AutoPas class as extern template instantiation. It is instantiated in AutoPasClass.cpp.
extern template class autopas::AutoPas<Particle>;

/**
 * Floating point precision for command line output.
 */
constexpr int floatStringPrecision = 3;

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
std::string timerToString(const std::string &name, long timeNS, int numberWidth = 0ul, long maxTime = 0ul) {
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

int main(int argc, char **argv) {
  // might be moved somewhere else if we add proper object orientation for the main logic
  struct Timers {
    autopas::utils::Timer total{};
    autopas::utils::Timer initialization{};
    autopas::utils::Timer simulation{};
    autopas::utils::Timer integrator{};
    autopas::utils::Timer collisionDetection{};
    autopas::utils::Timer containerUpdate{};
    autopas::utils::Timer output{};
  } __attribute__((aligned(128))) timers;

  timers.total.start();
  timers.initialization.start();

  Logger logger;

  logger.get()->set_level(spdlog::level::debug);

  // Default config path
  const auto *cfgFilePath = LoadConfig::defaultCfgPath;

  // Read in config if given
  if (argc > 1) cfgFilePath = argv[1];
  const auto config = LoadConfig::loadConfig(cfgFilePath, logger);
  logger.log(Logger::Level::info, "Config loaded.");

  const auto iterations = config["sim"]["iterations"].as<size_t>();
  const auto vtkWriteFrequency = config["io"]["vtkWriteFrequency"].as<size_t>();

  using AutoPas_t = autopas::AutoPas<Particle>;

  // initialization of autopas
  AutoPas_t autopas;
  const auto maxAltitude = config["sim"]["maxAltitude"].as<double>();
  const auto cutoff = config["autopas"]["cutoff"].as<double>();
  const auto desiredCellsPerDimension = config["autopas"]["desiredCellsPerDimension"].as<double>();

  autopas.setBoxMin({-maxAltitude, -maxAltitude, -maxAltitude});
  autopas.setBoxMax({maxAltitude, maxAltitude, maxAltitude});
  autopas.setCutoff(cutoff);
  // Set the size (relative to cutoff) of the cells so that roughly the desired number of cells per dimension is reached
  autopas.setCellSizeFactor((maxAltitude * 2.) / (cutoff * desiredCellsPerDimension));
  autopas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autopas.setAllowedContainers({autopas::ContainerOption::verletListsCells});
  autopas.init();

  // initialization of the integrator
  std::array<bool, 8> selectedPropagatorComponents{config["sim"]["prop"]["useKEPComponent"].as<bool>(),
                                                   config["sim"]["prop"]["useJ2Component"].as<bool>(),
                                                   config["sim"]["prop"]["useC22Component"].as<bool>(),
                                                   config["sim"]["prop"]["useS22Component"].as<bool>(),
                                                   config["sim"]["prop"]["useSOLComponent"].as<bool>(),
                                                   config["sim"]["prop"]["useLUNComponent"].as<bool>(),
                                                   config["sim"]["prop"]["useSRPComponent"].as<bool>(),
                                                   config["sim"]["prop"]["useDRAGComponent"].as<bool>()};

  auto fo = std::make_shared<FileOutput<AutoPas_t, Particle>>(
      autopas, config["io"]["output_file"].as<std::string>(), OutputFile::CSV, selectedPropagatorComponents);
  auto accumulator = std::make_shared<Acceleration::AccelerationAccumulator<AutoPas_t, Particle>>(
      selectedPropagatorComponents, autopas, 0.0, *fo);
  auto integrator = std::make_shared<Integrator<AutoPas_t, Particle>>(autopas, *accumulator, 1e-1);

  // Read in scenario
  auto actualSatellites =
      DatasetReader::readDataset(std::string(DATADIR) + config["io"]["posFileName"].as<std::string>(),
                                 std::string(DATADIR) + config["io"]["velFileName"].as<std::string>());
  logger.log(Logger::Level::debug, "Parsed {} satellites", actualSatellites.size());

  double minAltitudeFound{std::numeric_limits<double>::max()};
  double maxAltitudeFound{0.};
  // Convert satellites to particles
  for (const auto &particle : actualSatellites) {
    auto pos = particle.getPosition();
    double altitude = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
    minAltitudeFound = std::max(minAltitudeFound, altitude);
    maxAltitudeFound = std::max(maxAltitudeFound, altitude);
    if (altitude < maxAltitude) {
      autopas.addParticle(particle);
    }
  }
  logger.log(Logger::Level::info, "Min altitude is {}", minAltitudeFound);
  logger.log(Logger::Level::info, "Max altitude is {}", maxAltitudeFound);

  logger.log(Logger::Level::info, "Number of particles: {}", autopas.getNumberOfParticles());

  timers.initialization.stop();

  // main-loop skeleton
  timers.simulation.start();
  for (size_t i = 0; i < iterations; ++i) {
    // update positions
    timers.integrator.start();
    integrator->integrate(false);
    timers.integrator.stop();

    // TODO MPI: handle particle exchange between ranks
    timers.containerUpdate.start();
    // (potentially) update the internal data structure and collect particles which are leaving the container.
    const auto [escapedParticles, containerUpdated] = autopas.updateContainer();
    timers.containerUpdate.stop();

    if (not escapedParticles.empty()) {
      logger.log(Logger::Level::err, "Particles are escaping! \n{}", escapedParticles);
    }
    // TODO Check for particles that burn up

    timers.collisionDetection.start();
    // pairwise interaction
    CollisionFunctor collisionFunctor(cutoff);
    autopas.iteratePairwise(&collisionFunctor);
    auto collisions = collisionFunctor.getCollisions();
    logger.log(Logger::Level::info, "Iteration {} - Close encounters: {}", i, collisions.size());
    for (const auto &[p1, p2] : collisions) {
      logger.log(Logger::Level::debug, "{} | {}", p1->getID(), p2->getID());
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
      vtkWriter.printResult(allParticles);
    }
    timers.output.stop();
  }
  timers.simulation.stop();
  timers.total.stop();

  // dump all timers
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
      "  Output                ", timers.output.getTotalTime(), maximumNumberOfDigits, timeTotal);
  std::cout << timerToString("One iteration           ", timeSim / iterations, maximumNumberOfDigits, timeTotal);

  return 0;
}
