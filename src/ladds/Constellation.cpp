//
// Created by albert on 12.11.21.
//

#include "Constellation.h"

#include <cstdlib>
#include <iostream>
#include <string>
Constellation::Constellation(const std::string &constellation_data_str, int interval) {
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
  // convert vector to deque (should be changed some time)
  constellationSize = static_cast<int>(sats.size());
  for (int i = 0; i < constellationSize; i++) {
    satellites.push_back(sats.at(i));
  }

  startTime = std::stoi(argStartTime);
  duration = std::stoi(argDuration);
  timeActive = 0;
  this->interval = interval;
  simulationTime = 0;

  status = Status::inactive;

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
    shells.push_back({altitude, inclination, nPlanes, satsPerPlane});
    std::getline(shellParameters, tmp_string);
  }
  shellParameters.close();

  // determine times when each shell has its deployment started
  double timestamp = 0;
  timestamps.push_back(0);
  for (auto shell : shells) {
    timestamp += (shell[2] * shell[3] / constellationSize) * duration;
    timestamps.push_back(timestamp);
  }

  // for each shell determine the interval a new plane is added, each shell has its own timestep
  for (int i = 0; i < timestamps.size() - 1; i++) {
    timeSteps.push_back((timestamps.at(i + 1) - timestamps.at(i)) / shells.at(i)[2]);  // = duration_i / nPlanes_i
  }

  currentShellIndex = 0;
  planesDeployed = 0;
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

      while (timeActive >= timestamps.at(currentShellIndex) + planesDeployed * timeSteps.at(currentShellIndex)) {
        int planeSize = static_cast<int>(shells.at(currentShellIndex)[3]);
        for (int i = 0; i < planeSize; i++) {
          particles.push_back(satellites.at(0));
          satellites.pop_front();
        }
        planesDeployed++;

        // if all planes of the shell are done increment currentShell and set planesDeployed to 0
        if (planesDeployed >= shells.at(currentShellIndex)[2]) {
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

int Constellation::getConstellationSize() const {
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