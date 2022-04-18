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

#include "ladds/io/ConfigReader.h"
#include "ladds/particle/Particle.h"

/**
 * The Constellation class contains a collection of Particles that inserts these particles into
 * the simulation over time based on .csv files and a .yaml file as created by the
 * ConstellationGeneration notebook
 */
class Constellation {
 public:
  /**
   * Constructs a constellation object
   * @param constellationConfig : ConfigReader object with the constellation data. Valid
   * constellation data can be created using the ConstellationGeneration notebook
   * and must be in the projects data folder
   * @param interval : the interval of satellites being added to the simulation is
   * passed for internal logic
   * @param altitudeDeviation : used to create satellites with normally distributed
   * altitudes. Equals the standard deviation of a normal distribution
   * @param coefficientOfDrag c_D used to initialize all satellites.
   */
  Constellation(ConfigReader &constellationConfig, ConfigReader &config);

    /**
    * sets internal attribute startTime according to the passed date string
    * startTime_str
    * @param startTime a point in time either in iterations or as a date string. if
    * the string represents a natural number, it is considered as an iteration
    * and the string is converted to a number, if it is a date string it is converted
    * to an iteration timestamp before startTime is set to that value
    */
    void setStartTime(const std::string &startTime_str, const std::string &refTime_str);

    /**
     * sets internal attribute duration according to the passed string parameter
     * duration_str
     * @param duration_str represents the duration of deployment in either iterations
     * or days. the parameter is considered as a count of days when its last character
     * equals 'd' and an iteration count otherwise
     */
    void setDuration(const std::string &duration_str);

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

  /**
   * getter for constellationName = name of this constellation
   * @return std::sting : constellationName
   */
  [[nodiscard]] std::string getConstellationName() const;

  /**
   * getter for startTime = timestamp (iteration) when constellation insertion starts
   * @return size_t : timestamp (iteration) when constellation insertion starts
   */
  [[nodiscard]] long getStartTime() const;

  /**
   * getter for duration = timespan over which constellation insertion takes place
   * @return size_t : timespan over which constellation insertion takes place
   */
  [[nodiscard]] size_t getDuration() const;

 private:
  /**
   * stores the satellites of the constellation that have not been added to the simulation
   */
  std::deque<Particle> satellites{};

  /**
   * the name of the constellation
   */
  std::string constellationName;

  /**
   * Reads the passed position and velocity csv files. Returns a vector of particles.
   * @param position_filepath
   * @param velocity_filepath
   * @param coefficientOfDrag
   * @return
   */
  static std::vector<Particle> readDatasetConstellation(const std::string &position_filepath,
                                                        const std::string &velocity_filepath,
                                                        double coefficientOfDrag);

  /**
   * changes the pos vector by adding a random, normal distributed offset to the altitude
   * (offset dependent on altitudeVariance)
   * @param pos input position
   * @return new position with random altitude
   */
  std::array<double, 3> randomDisplacement(const std::array<double, 3> &pos);

  /**
   * parses datestring with expected format <year>/<month>/<day>
   * @param datestr datestring
   * @return integer array with year, month, and day
   */
   static std::array<int,3> parseDatestring(const std::string &date_str);

  /**
   * iteration from which constellation starts being added to the simulation
   */
  long startTime = 0;

  /**
   * time span over which satellites of the constellation are being added
   */
  size_t duration = 0;

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
   * deltaT of the simulation passed to constellations for converting datestring time
   * into simulation time expressed in iterations
   */
  double deltaT;

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
