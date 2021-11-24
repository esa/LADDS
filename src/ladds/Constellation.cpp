//
// Created by albert on 12.11.21.
//

#include "Constellation.h"
#include <cstdlib>
#include <iostream>
#include <string>
Constellation::Constellation(const std::string &constellation, int interval) {
  //split the 3 comma seperated arguments
  int seperator1 = constellation.find(',',0);
  int seperator2 = constellation.find(',',seperator1+1);
  std::string arg1 = constellation.substr(0,seperator1);                                           //  arg1     s1  arg2  s2   arg3    // n: 4
  std::string arg2 = constellation.substr(seperator1+1,seperator2-seperator1-1);                //a  b  c  d  ,  e  f  ,  g  h  i   // n: 7-4-1=2
  std::string arg3 = constellation.substr(seperator2+1,constellation.size()-seperator2-1);      //0  1  2  3  4  5  6  7  8  9  10  // n: 11-7-1=3

  //set variables using 3 args
  std::string path = arg1;
  std::vector<Particle> sats = readDatasetConstellation(std::string(DATADIR) + path + "/pos.csv",std::string(DATADIR) + path + "/v.csv");
  //convert vector to deque (should be changed some time)
  constellationSize = sats.size();
  for(int i = 0;i<constellationSize;i++) {
    satellites.push_back(sats.at(i));
  }

  startTime = std::stoi(arg2);
  duration = std::stoi(arg3);
  timeActive = 0;
  this->interval = interval;
  simulationTime = 0;

  status = 'i';

  std::ifstream shellParameters(std::string(DATADIR) + arg1 + "/shells.txt");
  std::string tmp_string = "";
  std::getline(shellParameters,tmp_string);
  double altitude,inclination,nPlanes,satsPerPlane;
  while(tmp_string != ""){
    std::istringstream numStream(tmp_string);
    numStream >> altitude;
    numStream >> inclination;
    numStream >> nPlanes;
    numStream >> satsPerPlane;
    shells.push_back({altitude,inclination,nPlanes,satsPerPlane});
    std::getline(shellParameters,tmp_string);
  }
  shellParameters.close();

  //determine times when each shell has its deployment started
  double timestamp = 0;
  timestamps.push_back(0);
  for(auto shell:shells){
    timestamp += (shell[2]*shell[3]/constellationSize)*duration;
    timestamps.push_back(timestamp);
  }

  //for each shell determine the interval a new plane is added, each shell has its own timestep
  for(int i = 0;i < timestamps.size()-1;i++) {
    timeSteps.push_back((timestamps.at(i+1)-timestamps.at(i))/shells.at(i)[2]);   // = duration_i / nPlanes_i
  }

  currentShellIndex = 0;
  planesDeployed = 0;
}

std::vector<Particle> Constellation::tick() {
  std::vector<Particle> particles{};
  switch(status) {
    case 'd':
      // do nothing
      break;
    case 'i':
      // check time and activate if startTime is reached
      if (simulationTime >= startTime) {
        status = 'a';
      } else {
        break;
      }
    case 'a':

      while (timeActive >= timestamps.at(currentShellIndex) + planesDeployed * timeSteps.at(currentShellIndex)) {

        int planeSize = shells.at(currentShellIndex)[3];
        for(int i = 0;i<planeSize;i++) {
          particles.push_back(satellites.at(0));
          satellites.pop_front();
        }
        planesDeployed++;

        //if all planes of the shell are done increment currentShell and set planesDeployed to 0
        if(planesDeployed >= shells.at(currentShellIndex)[2]) {
          currentShellIndex++;
          //end the operation, if every shell has been deployed = set constellation to 'd' = deployed
          if(currentShellIndex >= shells.size()){
            status = 'd';
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




std::vector<Particle> Constellation::readDatasetConstellation(const std::string &position_filepath, const std::string &velocity_filepath) {
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