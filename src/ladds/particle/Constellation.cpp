/**
 * @file Constellation.cpp
 * @author albert
 * @date 12.11.21
 */

#include "Constellation.h"

#include <breakupModel/input/CSVReader.h>

#include <cstdlib>
#include <iostream>
#include <string>

int Constellation::seed = 10;

Constellation::Constellation(const YAML::Node &constellationConfig, size_t interval, double altitudeDeviation)
    : interval(interval), altitudeDeviation(altitudeDeviation) {
  std::string constellationName = constellationConfig["constellation"]["name"].as<std::string>();

  // set variables using 3 args
  std::vector<Particle> sats =
      readDatasetConstellation(std::string(DATADIR) + constellationName + "/pos_" + constellationName + ".csv",
                               std::string(DATADIR) + constellationName + "/v_" + constellationName + ".csv");

  generator.seed(seed++);
  distribution = std::normal_distribution<double>(0, this->altitudeDeviation);
  // convert vector to deque
  constellationSize = sats.size();
  for (size_t i = 0ul; i < constellationSize; ++i) {
    sats[i].setPosition(randomDisplacement(sats[i].getPosition()));
    satellites.push_back(sats[i]);
  }

  startTime = constellationConfig["constellation"]["startTime"].as<int>();
  duration = constellationConfig["constellation"]["duration"].as<int>();

  int nShells = constellationConfig["constellation"]["nShells"].as<int>();
  for (int i = 1; i <= nShells; i++) {
    std::string attribute = "shell" + std::to_string(i);
    shells.emplace_back<std::array<double, 4>>({constellationConfig[attribute]["altitude"].as<double>(),
                                                constellationConfig[attribute]["inclination"].as<double>(),
                                                constellationConfig[attribute]["nPlanes"].as<double>(),
                                                constellationConfig[attribute]["nSats"].as<double>()});
  }

  // determine times when each shell has its deployment started
  double timestamp = 0;
  timestamps.push_back(0);
  for (auto [alt, i, planes, nSats] : shells) {
    timestamp += (planes * nSats / static_cast<double>(constellationSize)) * duration;
    timestamps.push_back(timestamp);
  }

  // for each shell determine the interval a new plane is added, each shell has its own timestep
  for (size_t i = 0ul; i < timestamps.size() - 1; ++i) {
    timeSteps.push_back((timestamps[i + 1] - timestamps[i]) / shells[i][2]);  // = duration_i / nPlanes_i
  }
}

std::vector<Particle> Constellation::tick() {
  std::vector<Particle> particles{};
  switch (status) {
    case Status::deployed:
      // do nothing
      break;
    case Status::inactive:
      // check time and activate if startTime is reached
      if (simulationTime >= startTime) {
        status = Status::active;
      } else {
        break;
      }
    case Status::active:

      while (static_cast<double>(timeActive) >=
             timestamps[currentShellIndex] + planesDeployed * timeSteps[currentShellIndex]) {
        int planeSize = static_cast<int>(shells[currentShellIndex][3]);
        particles.reserve(planeSize);
        for (int i = 0; i < planeSize; i++) {
          particles.push_back(satellites[0]);
          satellites.pop_front();
        }
        planesDeployed++;

        // if all planes of the shell are done increment currentShell and set planesDeployed to 0
        if (planesDeployed >= shells[currentShellIndex][2]) {
          currentShellIndex++;
          // end the operation, if every shell has been deployed = set constellation to 'd' = deployed
          if (currentShellIndex >= shells.size()) {
            status = Status::deployed;
            break;
          }

          planesDeployed = 0;
        }
      }
      timeActive += interval;
      break;
  }
  simulationTime += interval;
  return particles;
}

size_t Constellation::getConstellationSize() const {
  return constellationSize;
}

std::vector<Particle> Constellation::readDatasetConstellation(const std::string &position_filepath,
                                                              const std::string &velocity_filepath) {
  CSVReader<double, double, double> pos_csvReader{position_filepath, false};
  CSVReader<double, double, double> vel_csvReader{velocity_filepath, false};
  std::vector<Particle> particleCollection;

  auto positions = pos_csvReader.getLines();
  auto velocities = vel_csvReader.getLines();

  if (positions.size() != velocities.size()) {
    std::cout << "Error: Position and velocity file have different number of lines." << std::endl;
    return particleCollection;
  }

  particleCollection.reserve(positions.size());

  size_t particleId = 0;
  std::transform(positions.begin(),
                 positions.end(),
                 velocities.begin(),
                 std::back_insert_iterator<std::vector<Particle>>(particleCollection),
                 [&](const auto &pos, const auto &vel) {
                   const auto &[x, y, z] = pos;
                   const auto &[vx, vy, vz] = vel;

                   const std::array<double, 3> posArray = {x, y, z};
                   const std::array<double, 3> velArray = {vx, vy, vz};
                   return Particle(posArray, velArray, particleId++);
                 });
  return particleCollection;
}

std::array<double, 3> Constellation::randomDisplacement(std::array<double, 3> pos) {
  // the position pos is already the same as the direction vector from origin to pos
  // u = 1/length*pos
  std::array<double, 3> unitVector = autopas::utils::ArrayMath::normalize(pos);
  //
  double offset = distribution(generator);
  // npos = pos + offset * u
  return autopas::utils::ArrayMath::add(pos, autopas::utils::ArrayMath::mulScalar(unitVector, offset));
}