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

size_t Constellation::idBase = 1000000;

std::mt19937 Constellation::generator{42};

Constellation::Constellation(ConfigReader &constellationConfig,ConfigReader &config):
      constellationName(constellationConfig.get<std::string>("constellation/name")),
      interval(config.get<size_t>("io/constellationFrequency",1)),
      deltaT(config.get<double>("sim/deltaT")),
      /* altitudeSpread = 3 * sigma */
      distribution(std::normal_distribution<double>(0,config.get<double>("io/altitudeSpread", 0.0) / 3.0))
{

  auto coefficientOfDrag = config.get<double>("sim/prop/coefficientOfDrag");

  //set constellations insertion start and duration
  setStartTime(constellationConfig.get<std::string>("constellation/startTime"),
               config.get<std::string>("sim/referenceTime","2022/01/01"));
  setDuration(constellationConfig.get<std::string>("constellation/duration"));

  std::vector<Particle> sats =
      readDatasetConstellation(std::string(DATADIR) + constellationName + "/pos_" + constellationName + ".csv",
                               std::string(DATADIR) + constellationName + "/v_" + constellationName + ".csv",
                               coefficientOfDrag);
  // convert vector to deque and add noise to altitude
  constellationSize = sats.size();
  for (size_t i = 0ul; i < constellationSize; ++i) {
    sats[i].setPosition(randomDisplacement(sats[i].getPosition()));
    satellites.push_back(sats[i]);
  }

  int nShells = constellationConfig.get<int>("constellation/nShells");
  for (int i = 1; i <= nShells; i++) {
    std::string attribute = "shell" + std::to_string(i);
    shells.emplace_back<std::array<double, 4>>({constellationConfig.get<double>(attribute + "/altitude"),
                                                constellationConfig.get<double>(attribute + "/inclination"),
                                                constellationConfig.get<double>(attribute + "/nPlanes"),
                                                constellationConfig.get<double>(attribute + "/nSats")});
  }

  // determine times when each shell has its deployment started
  double timestamp = 0;
  std::vector<double> timestamps;
  timestamps.reserve(nShells+1);
  timestamps.push_back(0);

  for (auto [alt, i, planes, nSats] : shells) {
    timestamp += (planes * nSats / static_cast<double>(constellationSize)) * static_cast<double>(duration);
    timestamps.push_back(timestamp);
  }
  std::cout << "schedule for " << constellationName << ":" << std::endl;
  schedule.resize(nShells);
  for(int i = 0ul; i < timestamps.size() - 1;++i){
    int nPlanes = static_cast<int>(shells[i][2]);
    double timeStepSize = (timestamps[i + 1] - timestamps[i]) / nPlanes;

    schedule[i].reserve(nPlanes);
    for(int j = 0ul;j < nPlanes;j++){
      schedule[i].push_back(timestamps[i]+j*timeStepSize);
      std::cout << timestamps[i]+j*timeStepSize << " ";
    }
    std::cout << std::endl;
  }

  idBase += 1000000;
}

void Constellation::setStartTime(const std::string &startTime_str, const std::string &refTime_str) {
    //date string
    if(startTime_str.find('/') != std::string::npos) {
        std::array<int,3> dateArrayStart = parseDatestring(startTime_str);
        std::array<int,3> dateArrayRef = parseDatestring(refTime_str);
        struct tm stm = {0,0,0,dateArrayStart[2],dateArrayStart[1] - 1,dateArrayStart[0]};
        struct tm t0 = {0,0,0,dateArrayRef[2],dateArrayRef[1] - 1,dateArrayRef[0]};

        time_t stime = std::mktime(&stm) - std::mktime(&t0);
        startTime = static_cast<long>(static_cast<double>(stime) / deltaT);
        return;
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

std::vector<Particle> Constellation::tick() {
  std::vector<Particle> particles{};
  switch (status) {
    case Status::deployed:
      // do nothing
      break;
    case Status::inactive:
      // check time and activate if startTime is reached
      if (static_cast<long>(simulationTime) >= startTime) {
        status = Status::active;
        timeActive = simulationTime - startTime;
        std::cout << "timeActive for " << constellationName << ": " << timeActive << std::endl;
      } else {
        break;
      }
    case Status::active:

      while (static_cast<double>(timeActive) >= schedule[currentShellIndex][planesDeployed]) {
        std::cout << "insertion  " << constellationName << ": " << timeActive << std::endl;
        auto planeSize = static_cast<size_t>(shells[currentShellIndex][3]);
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
  simulationTime += interval;
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

std::vector<Particle> Constellation::readDatasetConstellation(const std::string &position_filepath,
                                                              const std::string &velocity_filepath,
                                                              double coefficientOfDrag) {
  size_t particleId = idBase;

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
                   const double mass{1.};
                   const double radius{1.};
                   // TODO: get proper constellation name
                   return Particle(posArray,
                                   velArray,
                                   particleId++,
                                   "Constellation",
                                   Particle::ActivityState::evasivePreserving,
                                   mass,
                                   radius,
                                   Particle::calculateBcInv(0., mass, radius, coefficientOfDrag));
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

std::array<int,3> Constellation::parseDatestring(const std::string &date_str) {
    std::array<int,3> dateArray{};
    std::string current_str = date_str;
    size_t current_idx = current_str.find('/');
    for(int i = 0;i<3;i++) {
        dateArray[i] = std::stoi(current_str.substr(0,current_idx));
        if(i != 2) {
            current_str = current_str.erase(0,current_idx + 1);
            current_idx = current_str.find('/');
        }
    }
    return dateArray;
}
