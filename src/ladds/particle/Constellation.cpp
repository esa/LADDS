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
#include <time.h>

size_t Constellation::particleID = 1000000;

std::mt19937 Constellation::generator{42};

Constellation::Constellation(const YAML::Node &constellationConfig,ConfigReader &config) {
  interval = config.get<int>("io/constellationFrequency",1,true);
  // altitudeSpread = 3 * sigma , altitudeDeviation = sigma (= standardDeviation)
  altitudeDeviation = config.get<double>("io/altitudeSpread") / 3.0;
  deltaT = config.get<double>("sim/deltaT");
  setStartTime(constellationConfig["constellation"]["startTime"].as<std::string>(),
          config.get<std::string>("sim/referenceTime","2022/01/01"));
  setDuration(constellationConfig["constellation"]["duration"].as<std::string>());

  constellationName = constellationConfig["constellation"]["name"].as<std::string>();

  // set variables using 3 args
  std::vector<Particle> sats =
      readDatasetConstellation(std::string(DATADIR) + constellationName + "/pos_" + constellationName + ".csv",
                               std::string(DATADIR) + constellationName + "/v_" + constellationName + ".csv");

  distribution = std::normal_distribution<double>(0, this->altitudeDeviation);
  // convert vector to deque
  constellationSize = sats.size();
  for (size_t i = 0ul; i < constellationSize; ++i) {
    sats[i].setPosition(randomDisplacement(sats[i].getPosition()));
    satellites.push_back(sats[i]);
  }

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
    timestamp += (planes * nSats / static_cast<double>(constellationSize)) * static_cast<double>(duration);
    timestamps.push_back(timestamp);
  }

  schedule.resize(nShells);
  for(int i = 0ul; i < timestamps.size() - 1;++i){
      int nPlanes = shells[i][2];
      double timeStepSize = (timestamps[i + 1] - timestamps[i]) / nPlanes;

      schedule[i].reserve(nPlanes);
      for(int j = 0ul;j < nPlanes;j++){
          schedule[i].push_back(timestamps[i]+j*timeStepSize);
      }
  }

  // prepare next ID base for next constellation (C1 starts at 1M, C2 starts at 2M ...)
  particleID = particleID + 1000000 - constellationSize;

  //when the simulation is loaded from a checkpoint the past satellites should NOT be added
  //(but satellites that start before 0 SHOULD be added when the negative startTime does not stem from a checkpoint)
  bool checkpointed = config.defines("io/checkpoint/file",false) &&
          (config.get<std::string>("io/checkpoint/file","") != "");

  if(checkpointed && startTime < 0) {
      status = Status::active;
      timeActive = -startTime;
      //simulate tick without adding anything until present is reached
      while(static_cast<double>(timeActive) > schedule[currentShellIndex][planesDeployed]){
          planesDeployed++;
          int planeSize = static_cast<int>(shells[currentShellIndex][3]);
          for (int i = 0; i < planeSize; i++) {
              satellites.pop_front();
          }
          if(planesDeployed >= shells[currentShellIndex][2]){
              currentShellIndex++;
              planesDeployed = 0;

              if(currentShellIndex > shells.size()){
                  status = Status::deployed;
                  break;
              }
          }
      }
  }
}

std::vector<Particle> Constellation::tick(size_t simulationTime) {
  std::vector<Particle> particles{};
  switch (status) {
    case Status::deployed:
      // do nothing
      break;
    case Status::inactive:
      // check time and activate if startTime is reached
      if (simulationTime >= startTime || startTime < 0) {
        status = Status::active;
        //if constellation has been scheduled before simulationStart, timeActive is set
        //accordingly to insert as much as is due
        timeActive = simulationTime - startTime;
      } else {
        break;
      }
    case Status::active:

      while (static_cast<double>(timeActive) >= schedule[currentShellIndex][planesDeployed]) {
        std::cout << constellationName << " " << schedule[currentShellIndex][planesDeployed] << std::endl;
        int planeSize = static_cast<int>(shells[currentShellIndex][3]);
        particles.reserve(particles.capacity()+planeSize);
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
  return particles;
}

size_t Constellation::getConstellationSize() const {
  return constellationSize;
}

std::string Constellation::getConstellationName() const {
    return constellationName;
}

long Constellation::getStartTime() const {
    return startTime;
}

size_t Constellation::getDuration() const {
    return duration;
}

void Constellation::setStartTime(const std::string &startTime_str, const std::string &referenceTime_str) {
    //date string
    if(startTime_str.size() > 4){
        if(startTime_str[4] == '/'){
            int year = std::stoi(startTime_str.substr(0, 4));
            int month = std::stoi(startTime_str.substr(5, 2)) - 1;
            int day = std::stoi(startTime_str.substr(8, 2));
            int refYear = std::stoi(referenceTime_str.substr(0, 4));
            int refMonth = std::stoi(referenceTime_str.substr(5, 2)) - 1;
            int refDay = std::stoi(referenceTime_str.substr(8, 2));
            struct tm stm = {0,0,0,day,month,year};
            struct tm t0 = {0,0,0,refDay,refMonth,refYear};

            time_t stime = std::mktime(&stm) - std::mktime(&t0);
            //integer division cutting off anything smaller than 1ms
            startTime = static_cast<long>(static_cast<long>(stime)*1000) / static_cast<long>(deltaT*1000.0);
            return;
        }
    }
    //iteration
    startTime = std::stoi(startTime_str);
}

void Constellation::setDuration(const std::string &duration_str) {
    if(duration_str[duration_str.size()-1] == 'd') {
        duration = static_cast<size_t>(24*60*60*std::stoi(duration_str.substr(0,duration_str.size()-1)) / deltaT);
    } else {
        duration = std::stoi(duration_str);
    }
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

  std::transform(positions.begin(),
                 positions.end(),
                 velocities.begin(),
                 std::back_insert_iterator<std::vector<Particle>>(particleCollection),
                 [&](const auto &pos, const auto &vel) {
                   const auto &[x, y, z] = pos;
                   const auto &[vx, vy, vz] = vel;

                   const std::array<double, 3> posArray = {x, y, z};
                   const std::array<double, 3> velArray = {vx, vy, vz};
                   return Particle(posArray, velArray, particleID++);
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