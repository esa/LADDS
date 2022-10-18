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

#include "ladds/distributedMemParallelization/DomainDecomposition.h"
#include "ladds/io/ConfigReader.h"
#include "ladds/particle/Particle.h"

namespace LADDS {

/**
 * The Constellation class contains a collection of Particles that inserts these into
 * the simulation over time, based on .csv files and a .yaml file as created by the
 * ConstellationGeneration notebook.
 */
class Constellation {
 public:
  /**
   * Constructs a constellation object
   * @param constellationConfig : ConfigReader object with the constellation data. Valid
   * constellation data can be created using the ConstellationGeneration notebook
   * and must be placed in the projects data folder.
   * @param config : ConfigReader object with the simulation information, reads constellationFrequency,
   * deltaT,coefficientOfDrag.
   */
  Constellation(ConfigReader &constellationConfig, ConfigReader &config);

  /**
   * Determines which satellites are being added to the simulation by adding each shell
   * within a time span proportional to the shells size. shells are added plane by plane
   * and linearly over time.
   * @param domainDecomposition
   * @return std::vector<Particle> : satellites to be added to the simulation.
   */
  std::vector<Particle> tick(DomainDecomposition &domainDecomposition);

  /**
   * Offsets all local constellation IDs by the parameter baseId to create global IDs.
   * (globalId = localId + baseId)
   * @param baseId
   */
  void moveConstellationIDs(const size_t baseId);

  /**
   * Getter for constellationSize = number of satellites in constellation.
   * @return int : constellationSize.
   */
  [[nodiscard]] size_t getConstellationSize() const;

  /**
   * Getter for constellationName = name of this constellation.
   * @return std::sting : constellationName.
   */
  [[nodiscard]] std::string getConstellationName() const;

  /**
   * Getter for startTime = timestamp (iteration) when constellation insertion starts.
   * @return size_t : timestamp (iteration) when constellation insertion starts.
   */
  [[nodiscard]] long getStartTime() const;

  /**
   * Getter for duration = timespan over which constellation insertion takes place.
   * @return size_t : timespan over which constellation insertion takes place.
   */
  [[nodiscard]] size_t getDuration() const;

 private:
  /**
   * Stores the satellites of the constellation that have not been added to the simulation.
   */
  std::deque<Particle> satellites{};

  /**
   * The name of the constellation.
   */
  std::string constellationName;

  /**
   * Reads the passed position and velocity csv files. Returns a vector of particles.
   * @param positionFilepath : path to csv file with the satellites' positions.
   * @param velocityFilepath : path to csv file with the satellites' velocities.
   * @param coefficientOfDrag : coefficient of drag.
   * @return vector with Particle objects.
   */
  static std::vector<Particle> readDatasetConstellation(const std::string &positionFilepath,
                                                        const std::string &velocityFilepath,
                                                        double coefficientOfDrag);

  /**
   * Sets internal attribute startTime according to the passed date string
   * startTimeStr.
   * @param startTime a point in time either in iterations or as a date string. if
   * the string represents a natural number, it is considered as an iteration
   * and the string is converted to a number, if it is a date string it is converted
   * to an iteration timestamp before startTime is set to that value.
   */
  void setStartTime(const std::string &startTimeStr, const std::string &refTimeStr);

  /**
   * Sets internal attribute duration according to the passed string parameter
   * durationStr.
   * @param durationStr represents the duration of deployment in either iterations
   * or days. the parameter is considered as a count of days when its last character
   * equals 'd' and an iteration count otherwise.
   */
  void setDuration(const std::string &durationStr);

  /**
   * Changes the pos vector by adding a random, normal distributed offset to the altitude
   * (offset dependent on altitudeVariance).
   * @param pos input position.
   * @return New position with random altitude.
   */
  std::array<double, 3> randomDisplacement(const std::array<double, 3> &pos);

  /**
   * Parses datestring with expected format <year>-<month>-<day>
   * @param datestr datestring.
   * @return Integer array with year, month, and day.
   */
  static std::array<int, 3> parseDatestring(const std::string &dateStr);

  /**
   * Iteration from which constellation starts being added to the simulation.
   */
  long startTime = 0;

  /**
   * Time span over which satellites of the constellation are being added. (iterations)
   */
  size_t duration = 0;

  /**
   * Internal clock that determines which satellites are added to the simulation,
   * starts the count when constellation is set to 'a' = active. (iterations)
   */
  size_t timeActive = 0;

  /**
   * The interval of satellites being added to the simulation is
   * passed for internal logic. (iterations)
   */
  size_t interval = 0;

  /**
   * deltaT of the simulation passed to constellations for converting datestring time
   * into simulation time expressed in iterations. (simulation seconds)
   */
  double deltaT;

  /**
   * The constellations state is set to 'a' = active whenever
   * simulationTime reaches startTime. (iterations)
   */
  size_t simulationTime = 0;

  /**
   * The three different possible internal states of a constellation object:
   * inactive: startTime has not been reached yet.
   * active: the constellation is currently being added to the simulation.
   * deployed: the constellation is fully deployed, and tick becomes a NOOP.
   */
  enum Status { inactive, active, deployed };

  /**
   * Variable that holds the internal state of the constellation object that determines
   * the behaviour of tick(). There are 3 different states:
   * inactive: startTime has not been reached yet.
   * active: the constellation is currently being added to the simulation.
   * deployed: the constellation is fully deployed, and tick becomes a NOOP.
   */
  Status status = Status::inactive;  // active , inactive ,deployed

  /**
   * Total size of the constellation. Does not change unlike satellites.size()
   */
  size_t constellationSize = 0ul;

  /**
   * Contains information of a shell: altitude, inclination, #planes, #satellitesPerPlane.
   */
  std::vector<std::array<double, 4>> shells{};

  /**
   * A vector containing each shells' schedule. if the entry at currentShellIndex and
   * planesDeployed is smaller than timeActive the corresponding plane is added.
   */
  std::vector<std::vector<double>> schedule{};

  /**
   * Keeps track of which shell will be added next.
   */
  size_t currentShellIndex = 0ul;

  /**
   * Keeps track of which plane will be added next.
   */
  int planesDeployed = 0;

  /**
   * Seeded/deterministic random number generator used to add noise to the
   * altitudes of satellites.
   */
  static std::mt19937 generator;

  /**
   * Normal distribution that determines the deviation of the satellites base
   * altitude. uses altitudeDeviation as parameter.
   */
  std::normal_distribution<double> distribution;
};

}  // namespace LADDS
