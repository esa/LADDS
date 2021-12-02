//
// Created by albert on 12.11.21.
//
#pragma once
#include <breakupModel/input/CSVReader.h>
#include <breakupModel/output/VTKWriter.h>

#include <array>
#include <deque>
#include <iostream>

#include "Particle.h"

/**
 * The Constellation class contains a collection of Particles that inserts these particles to
 * the simulation over time based on parameters a Constellation object is constructed with.
 * The object is constructed from a data string that consists of the comma seperated
 * arguments: directory path (directory with constellation information files), start time
 * (time when constellation is inserted), duration (time span of insertion)
 */
class Constellation {
 public:
  /**
   * Constructs a constellation object
   * @param constellation : string formatted as a comma seperated 3-tuple containing
   * a path to a directory with constellation information, the time the constellation
   * deployment is started, and the duration of the duration
   * @param interval : the interval of satellites being added to the simulation is
   * passed for internal logic
   */
  Constellation(const std::string &constellation_data_str, int interval);
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
  std::deque<Particle> satellites;

  /**
   * Reads the passed position and velocity csv files. Returns a vector of particles.
   */
  static std::vector<Particle> readDatasetConstellation(const std::string &position_filepath,
                                                        const std::string &velocity_filepath);

  /**
   * iteration from which constellation starts being added to the simulation
   */
  int startTime;

  /**
   * time span over which satellites of the constellation are being added
   */
  int duration;

  /**
   * internal clock that determines which satellites are added to the simulation,
   * starts the count when constellation is set to 'a' = active
   */
  int timeActive = 0;

  /**
   * the interval of satellites being added to the simulation is
   * passed for internal logic
   */
  int interval;

  /**
   * multiples of interval. the constellations state is set to 'a' = active whenever
   * simulationTime reaches startTime
   */
  int simulationTime = 0;

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
  size_t constellationSize;

  /**
   * contains information of a shell: altitude, inclination, #planes, #satellitesPerPlane
   */
  std::vector<std::array<double, 4>> shells;

  /**
   * contains the time shell i begins its deployment at vector index i
   */
  std::vector<double> timestamps;

  /**
   * contains the time steps of shell i to enable adding each plane of
   * shell i linearly over time at vector index i
   */
  std::vector<double> timeSteps;

  /**
   * keeps track of which shell will be added next
   */
  int currentShellIndex;

  /**
   * keeps track of which plane will be added next
   */
  int planesDeployed;
};
