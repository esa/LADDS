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
    setStartTime(constellationConfig["constellation"]["startTime"].as<std::string>());
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

    auto nShells = constellationConfig["constellation"]["nShells"].as<size_t>();
    for (int i = 1; i <= nShells; i++) {
        std::string attribute = "shell" + std::to_string(i);
        shells.emplace_back<std::array<double, 4>>({constellationConfig[attribute]["altitude"].as<double>(),
                                                    constellationConfig[attribute]["inclination"].as<double>(),
                                                    constellationConfig[attribute]["nPlanes"].as<double>(),
                                                    constellationConfig[attribute]["nSats"].as<double>()});
    }
    calculateRevolutionTimes(nShells);

    // determine times when each shell has its deployment started
    double timestamp = 0;
    timestamps.push_back(0);
    for (auto [alt, i, planes, nSats] : shells) {
        timestamp += (planes * nSats / static_cast<double>(constellationSize)) * static_cast<double>(duration);
        timestamps.push_back(timestamp);
    }

    // for each shell determine the interval a new plane is added, each shell has its own timestep
    for (size_t i = 0ul; i < timestamps.size() - 1; ++i) {
        timeSteps.push_back((timestamps[i + 1] - timestamps[i]) / shells[i][2]);  // = duration_i / nPlanes_i
    }

    // prepare next ID base for next constellation (C1 starts at 1M, C2 starts at 2M ...)
    particleID = particleID + 1000000 - constellationSize;
}

std::vector<Particle> Constellation::tick() {

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
            if(planesDeployed != 0) {
                phaseShift = fmod(phaseShift + static_cast<double>(interval),
                                  revolutionTimes[currentShellIndex] / deltaT);
            }
            while (static_cast<double>(timeActive) >=
                   timestamps[currentShellIndex] + planesDeployed * timeSteps[currentShellIndex]) {
                int planeSize = static_cast<int>(shells[currentShellIndex][3]);
                particles.reserve(planeSize);
                std::array<std::array<double,3>,3> rotationMatrix = calculateRotationMatrix(satellites[0].getPosition(),satellites[0].getVelocity(),currentShellIndex);
                for (int i = 0; i < planeSize; i++) {
                    Particle nParticle = transformSatellite(satellites[0],rotationMatrix);
                    particles.push_back(nParticle);
                    satellites.pop_front();
                }
                planesDeployed++;

                // if all planes of the shell are done increment currentShell and set planesDeployed to 0
                if (planesDeployed >= shells[currentShellIndex][2]) {
                    phaseShift = 0;
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

void Constellation::setStartTime(const std::string &startTime_str) {
    //date string
    if(startTime_str.size() > 4){
        if(startTime_str[4] == '/'){
            int year = std::stoi(startTime_str.substr(0, 4));
            int month = std::stoi(startTime_str.substr(5, 2)) - 1;
            int day = std::stoi(startTime_str.substr(8, 2));
            struct tm stm = {0,0,0,day,month,year};
            struct tm t0 = {0,0,0,1,0,2022};
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

void Constellation::calculateRevolutionTimes(size_t nShells) {
    size_t satelliteIndex = 0;
    for(size_t i = 0;i<nShells;i++){
        Particle p = satellites[satelliteIndex];
        //t = 2*pi*r/|v| (kms to travel / kms per second)
        double revolutionTime = doublePi*(shells[i][0]+6371.0)/autopas::utils::ArrayMath::L2Norm(p.getVelocity());
        //double revolutionTime = doublePi*sqrt(pow(shells[i][0]+6371.0,3)/3.986004407799724e+5);
        revolutionTimes.push_back(revolutionTime);
        //first satellite of next shell (index += nPlanes*nSats)
        satelliteIndex += shells[i][2]*shells[i][3];
    }
}

std::array<std::array<double,3>,3> Constellation::calculateRotationMatrix(std::array<double,3> pos, std::array<double,3> v, size_t shellIndex){
    //https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle

    //1 e_x
    std::array<double,3> v1 = autopas::utils::ArrayMath::normalize(pos);
    //2 e_y
    std::array<double,3> v2 = autopas::utils::ArrayMath::normalize(v);
    //3 e_z
    std::array<double,3> u = {v1[2]*v2[1]-v1[1]*v2[2],v1[0]*v2[2]-v1[2]*v2[0],v1[1]*v2[0]-v1[0]*v2[1]};

    //std::cout << "PhaseShift: " << phaseShift << std::endl;

    //rotation angle, based on the proportion of passed time w.r.t. one period
    double angle = - doublePi * phaseShift * deltaT/revolutionTimes[shellIndex];

    //std::cout << "Angle: " << angle << std::endl;

    std::array<double,3> row1 = {cos(angle) + u[0]*u[0]* (1 - cos(angle)),
                                 u[0]*u[1]*(1-cos(angle))- u[2]*sin(angle),
                                 u[0]*u[2]*(1-cos(angle)) + u[1]*sin(angle)};
    std::array<double,3> row2 = {u[1]*u[0]*(1-cos(angle)) + u[2]*sin(angle),
                                 cos(angle)+ u[1]*u[1]*(1-cos(angle)),
                                 u[1]*u[2]*(1-cos(angle)) - u[0]*sin(angle)};
    std::array<double,3> row3 = {u[2]*u[0]*(1-cos(angle)) - u[1]*sin(angle),
                                 u[2]*u[1]*(1-cos(angle)) + u[0]*sin(angle),
                                 cos(angle)+ u[2]*u[2]*(1-cos(angle))};

    std::array<std::array<double,3>,3> A = {row1,row2,row3};

    return A;
}

Particle Constellation::transformSatellite(Particle p,const std::array<std::array<double,3>,3>& rotation){
    std::array<double,3> pos = p.getPosition();
    std::array<double,3> v = p.getVelocity();
    //1 pos' = rotationMatrix*pos
    std::array<double,3> npos{};
    npos[0] = autopas::utils::ArrayMath::dot(rotation[0],pos);
    npos[1] = autopas::utils::ArrayMath::dot(rotation[1],pos);
    npos[2] = autopas::utils::ArrayMath::dot(rotation[2],pos);
    //2 v' = rotationMatrix*v
    std::array<double,3> nv{};
    nv[0] = autopas::utils::ArrayMath::dot(rotation[0],v);
    nv[1] = autopas::utils::ArrayMath::dot(rotation[1],v);
    nv[2] = autopas::utils::ArrayMath::dot(rotation[2],v);

    /*
    double angleCalc = acos(autopas::utils::ArrayMath::dot(autopas::utils::ArrayMath::normalize(pos),autopas::utils::ArrayMath::normalize(npos)));
    double angleCalc2 = acos(autopas::utils::ArrayMath::dot(autopas::utils::ArrayMath::normalize(npos),autopas::utils::ArrayMath::normalize(nv)));
    double velocityDiff = autopas::utils::ArrayMath::L2Norm(v)-autopas::utils::ArrayMath::L2Norm(nv);
    std::cout << "angle between pos - npos: " << angleCalc << std::endl;
    std::cout << "angle between npos - nv: " << angleCalc2 << std::endl;
    std::cout << "difference in velocity magnitude: " << velocityDiff << std::endl;
    */
    p.setPosition(npos);
    p.setVelocity(nv);
    return p;
}