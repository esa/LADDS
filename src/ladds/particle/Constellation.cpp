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

std::mt19937 Constellation::generator{42};

Constellation::Constellation(const std::string &constellation_data_str, size_t interval, double altitudeDeviation)
    : interval(interval), altitudeDeviation(altitudeDeviation) {
  // split the 3 comma seperated arguments
  auto seperator1 = constellation_data_str.find(',', 0);
  auto seperator2 = constellation_data_str.find(',', seperator1 + 1);

  // argDirPath s1 argStartTime s2  argDuration    // n: 4
  // a   b   c   d   ,   e   f   ,   g   h   i   // n: 7-4-1=2
  // 0   1   2   3   4   5   6   7   8   9   10  // n: 11-7-1=3
  std::string argConstellationName = constellation_data_str.substr(0, seperator1);
  std::string argStartTime = constellation_data_str.substr(seperator1 + 1, seperator2 - seperator1 - 1);
  std::string argDuration =
      constellation_data_str.substr(seperator2 + 1, constellation_data_str.size() - seperator2 - 1);

  // set variables using 3 args
  std::vector<Particle> sats =
      readDatasetConstellation(std::string(DATADIR) + argConstellationName + "/pos_" + argConstellationName + ".csv",
                               std::string(DATADIR) + argConstellationName + "/v_" + argConstellationName + ".csv");

  distribution = std::normal_distribution<double>(0, this->altitudeDeviation);
  // convert vector to deque
  constellationSize = sats.size();
  for (size_t i = 0ul; i < constellationSize; ++i) {
    sats[i].setPosition(randomDisplacement(sats[i].getPosition()));
    satellites.push_back(sats[i]);
  }

  startTime = std::stoi(argStartTime);
  duration = std::stoi(argDuration);

  std::ifstream shellParameters(std::string(DATADIR) + argConstellationName + "/shells_" + argConstellationName +
                                ".txt");
  std::string tmp_string;
  std::getline(shellParameters, tmp_string);
  double altitude, inclination, nPlanes, satsPerPlane;
  while (!tmp_string.empty()) {
    std::istringstream numStream(tmp_string);
    numStream >> altitude;
    numStream >> inclination;
    numStream >> nPlanes;
    numStream >> satsPerPlane;
    shells.emplace_back<std::array<double, 4>>({altitude, inclination, nPlanes, satsPerPlane});
    std::getline(shellParameters, tmp_string);
  }
  shellParameters.close();

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

std::array<double, 3> Constellation::randomDisplacement(const std::array<double, 3> &pos) {
  // the position pos is already the same as the direction vector from origin to pos
  // u = 1/length*pos
  std::array<double, 3> unitVector = autopas::utils::ArrayMath::normalize(pos);
  //
  double offset = distribution(generator);
  // npos = pos + offset * u
  return autopas::utils::ArrayMath::add(pos, autopas::utils::ArrayMath::mulScalar(unitVector, offset));
}