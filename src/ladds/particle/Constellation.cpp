/**
 * @file Constellation.cpp
 * @author albert
 * @date 12.11.21
 */

#include "Constellation.h"

#include <breakupModel/input/CSVReader.h>

#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace LADDS {

std::mt19937 Constellation::generator{42};

Constellation::Constellation(ConfigReader &constellationConfig, ConfigReader &config)
    : constellationName(constellationConfig.get<std::string>("constellation/name")),
      interval(config.get<size_t>("io/constellationFrequency", 1)),
      deltaT(config.get<double>("sim/deltaT")),
      /* altitudeSpread = 3 * sigma */
      distribution(std::normal_distribution<double>(0, config.get<double>("io/altitudeSpread", 0.0) / 3.0)) {
  const auto coefficientOfDrag = config.get<double>("sim/prop/coefficientOfDrag");

  // set constellations insertion start and duration
  setStartTime(constellationConfig.get<std::string>("constellation/startTime"),
               config.get<std::string>("sim/referenceTime"));
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

  const auto nShells = constellationConfig.get<int>("constellation/nShells");
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
  timestamps.reserve(nShells + 1);
  timestamps.push_back(0);

  for (const auto [alt, i, planes, nSats] : shells) {
    timestamp += (planes * nSats / static_cast<double>(constellationSize)) * static_cast<double>(duration);
    timestamps.push_back(timestamp);
  }

  // calculate schedule with launch times for each plane in constellation
  schedule.resize(nShells);
  for (size_t i = 0ul; i < timestamps.size() - 1; ++i) {
    int nPlanes = static_cast<int>(shells[i][2]);
    double timeStepSize = (timestamps[i + 1] - timestamps[i]) / nPlanes;

    schedule[i].reserve(nPlanes);
    for (int j = 0ul; j < nPlanes; j++) {
      schedule[i].push_back(timestamps[i] + j * timeStepSize);
    }
  }
}

void Constellation::setStartTime(const std::string &startTimeStr, const std::string &refTimeStr) {
  // date string
  if (startTimeStr.find('-') != std::string::npos) {
    std::array<int, 3> dateArrayStart = parseDatestring(startTimeStr);
    std::array<int, 3> dateArrayRef = parseDatestring(refTimeStr);
    struct tm stm = {0, 0, 0, dateArrayStart[2], dateArrayStart[1] - 1, dateArrayStart[0]};
    struct tm t0 = {0, 0, 0, dateArrayRef[2], dateArrayRef[1] - 1, dateArrayRef[0]};

    time_t stime = std::mktime(&stm) - std::mktime(&t0);
    startTime = static_cast<long>(static_cast<double>(stime) / deltaT);
    return;
  }
  // iteration
  startTime = std::stoi(startTimeStr);
}

void Constellation::setDuration(const std::string &durationStr) {
  if (durationStr[durationStr.size() - 1] == 'd') {
    duration = static_cast<size_t>(24 * 60 * 60 * std::stoi(durationStr.substr(0, durationStr.size() - 1)) / deltaT);
  } else {
    duration = std::stoi(durationStr);
  }
}

std::vector<Particle> Constellation::tick(DomainDecomposition &domainDecomposition) {
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
      } else {
        break;
      }
    case Status::active:
      const auto [shellAltitude, shellInclination, planeNumber, planeSize] = shells[currentShellIndex];
      int myRank{};
      autopas::AutoPas_MPI_Comm_rank(domainDecomposition.getCommunicator(), &myRank);
      // inserting scheduled particles
      while (static_cast<double>(timeActive) >= schedule[currentShellIndex][planesDeployed]) {
        particles.reserve(particles.capacity() + static_cast<size_t>(planeSize));
        for (size_t i = 0; i < planeSize; i++) {
          // make sure only satellites that are on this rank are added.
          if (domainDecomposition.getRank(satellites.front().getPosition()) == myRank) {
            particles.push_back(satellites.front());
          }
          satellites.pop_front();
        }
        planesDeployed++;

        // if all planes of the shell are done increment currentShell and set planesDeployed to 0
        if (planesDeployed >= planeNumber) {
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

void Constellation::moveConstellationIDs(const size_t baseId) {
  for (auto &satellite : satellites) {
    satellite.setID(baseId + satellite.getID());
  }
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
  size_t particleId = 0;

  CSVReader<double, double, double> pos_csvReader{position_filepath, false};
  CSVReader<double, double, double> vel_csvReader{velocity_filepath, false};
  std::vector<Particle> particleCollection;

  auto positions = pos_csvReader.getLines();
  auto velocities = vel_csvReader.getLines();

  if (positions.size() != velocities.size()) {
    throw std::runtime_error("Error: Position and velocity file have different number of lines.");
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
                                   Particle::calculateBcInv(0., mass, radius, coefficientOfDrag),
                                   std::numeric_limits<size_t>::max());
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

std::array<int, 3> Constellation::parseDatestring(const std::string &dateStr) {
  std::array<int, 3> dateArray{};
  std::string currentStr = dateStr;
  for (auto &dateField : dateArray) {
    const size_t currentIdx = std::min(currentStr.find('-'), currentStr.size());
    dateField = std::stoi(currentStr.substr(0, currentIdx));
    currentStr = currentStr.erase(0, currentIdx + 1);
  }
  return dateArray;
}

}  // namespace LADDS