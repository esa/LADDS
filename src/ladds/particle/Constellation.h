/**
 * @file Constellation.h
 * @author albert
 * @date 12.11.21
 */
#pragma once

#include <yaml-cpp/yaml.h>

#include <array>
#include <deque>
#include <iostream>
#include <random>

#include "ladds/particle/Particle.h"

/**
 * The Constellation class contains a collection of Particles that inserts these particles into
 * the simulation over time based on .csv files and a .yaml file as created by the
 * ConstellationGeneration notebook +
 */
class Constellation {
 public:
  /**
   * Constructs a constellation object

   * @param constellationConfig : YAML::Node object with the constellation data. Valid
   * constellation data can be created using the ConstellationGeneration notebook
   * and must be in the projects data folder
   * @param interval : the interval of satellites being added to the simulation is
   * passed for internal logic
   * @param altitudeDeviation : used to create satellites with normally distributed
   * altitudes. Equals the standard deviation of a normal distribution
   */
  Constellation(const YAML::Node &constellationConfig, size_t interval, double altitudeDeviation);
  /**
   * determines which satellites are being added to the simulation by adding each shell
   * within a time span proportional to the shells size. shells are added plane by plane
   * and linearly over time
   * @return std::vector<Particle> : satellites to be added to the simulation
   */
  std::vector<Particle> tick();

  /**
   * getter for constellationSize = number of satellites in constellation
   * @return int : constellationSize
   */
  [[nodiscard]] size_t getConstellationSize() const;

 private:
  /**
   * stores the satellites of the constellation that have not been added to the simulation
   */
  std::deque<Particle> satellites{};

  /**
   * Reads the passed position and velocity csv files. Returns a vector of particles.
   */
  static std::vector<Particle> readDatasetConstellation(const std::string &position_filepath,
                                                        const std::string &velocity_filepath);

  /**
   * changes the pos vector by adding a random, normal distributed offset to the altitude
   * (offset dependent on altitudeVariance)
   * @param pos input position
   * @return new position with random altitude
   */
  std::array<double, 3> randomDisplacement(const std::array<double, 3> &pos);

  /**
   * iteration from which constellation starts being added to the simulation
   */
  int startTime = 0;

  /**
   * time span over which satellites of the constellation are being added
   */
  int duration = 0;

  /**
   * internal clock that determines which satellites are added to the simulation,
   * starts the count when constellation is set to 'a' = active
   */
  size_t timeActive = 0;

  /**
   * the interval of satellites being added to the simulation is
   * passed for internal logic
   */
  size_t interval = 0;

  /**
   * multiples of interval. the constellations state is set to 'a' = active whenever
   * simulationTime reaches startTime
   */
  size_t simulationTime = 0;

  /**
   * The three different possible internal states of a constellation object:
   * inactive: startTime has not been reached yet
   * active: the constellation is currently being added to the simulation
   * deployed: the constellation is fully deployed, and tick becomes a NOOP
   */
  enum Status { inactive, active, deployed };

  /**
   * variable that holds the internal state of the constellation object that determines
   * the behaviour of tick(). There are 3 different states:
   * inactive: startTime has not been reached yet
   * active: the constellation is currently being added to the simulation
   * deployed: the constellation is fully deployed, and tick becomes a NOOP
   */
  Status status = Status::inactive;  // active , inactive ,deployed

  /**
   * size of the constellation for internal use
   */
  size_t constellationSize = 0ul;

  /**
   * contains information of a shell: altitude, inclination, #planes, #satellitesPerPlane
   */
  std::vector<std::array<double, 4>> shells{};

  /**
   * contains the time shell i begins its deployment at vector index i
   */
  std::vector<double> timestamps{};

  /**
   * contains the time steps of shell i to enable adding each plane of
   * shell i linearly over time at vector index i
   */
  std::vector<double> timeSteps{};

  /**
   * keeps track of which shell will be added next
   */
  size_t currentShellIndex = 0ul;

  /**
   * keeps track of which plane will be added next
   */
  int planesDeployed = 0;

  /**
   * deviation parameter of the normal distribution that determines the deviation
   * of the satellites base altitude. Equals the standard deviation of a normal distribution
   */
  double altitudeDeviation;

  /**
   * seeded/deterministic random number generator used to add noise to the
   * altitudes of satellites
   */
  static std::mt19937 generator;

  /**
   * normal distribution that determines the deviation of the satellites base
   * altitude. uses altitudeDeviation as parameter
   */
  std::normal_distribution<double> distribution;
};
