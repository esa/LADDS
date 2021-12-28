# LADDS - Large-scale Deterministic Debris Simulation

Codebase for the ARIADNA Study between TU Munich and ESA's Advanced Concepts Team. A more detailed project description can be found on the [Advanced Concept Team's webpage](https://www.esa.int/gsp/ACT/projects/debris_hpc/).

## Requirements
* CMake >= 3.14
* make (build-essentials, or equivalent)
* A C++17 Compiler (recommended: gcc >=7 or clang >=8, only clang-10 is tested)
* OpenMP >= 4.5
* git (for fetching dependencies)
* [TBB](https://github.com/oneapi-src/oneTBB) (Breakup-Model needs this)

### Optional
* Doxygen
* clang-format-9

## Building
```bash
mkdir build && cd build
CC=clang CXX=clang++ ccmake ..  # Set Variables according to your preferences
make -j12                       # choose number according to your CPU
```

## Testing
Testing is done with help of [GoogleTest](https://github.com/google/googletest), which is downloaded by CMake.
```bash
cmake -DLADDS_BUILD_TESTS=ON .. # Should be enabled by default
make ladds_tests -j12
ctest -j12
```

## Calibrating AutoPas
### Enable and Analyze Tuning
If you are unsure what algorithmic configuration you want to use for AutoPas just let AutoPas guide you.
For this, two things need to be activated:
* In the `yaml` file:
```yaml
autopas:
  tuningMode: true
```
* The CMake variables:  `AUTOPAS_LOG_TUNINGDATA=ON` and `AUTOPAS_LOG_TUNINGRESULTS=ON`.

This will result in AutoPas testing all reasonable configurations and dumping the results in two `csv` files.
* `AutoPas_tuningData.csv` contains the timing data of all samples AutoPas collected.
* `AutoPas_tuningResults.csv` contains the result of each tuning phase.
  If this file is empty, you either forgot to activate tuning mode or did not run enough iterations for a 
  tuning cycle to finish. Currently, about `100*rebuildFrequency` iterations are needed.

### Configuring AutoPas manually
When `tuningMode` is disabled it is possible to manually select the Algorithms AutoPas should use through
the `yaml` file. See `default_cfg.yaml` for the syntax.

## Processing TLE Input 
Data on current satellites etc. is often found [online](https://www.space-track.org/) in the [TLE format](https://en.wikipedia.org/wiki/Two-line_element_set). We include a Jupyter notebook which can be used to process TLE data with pykep to create and analyze suitable datasets. Detailed instructions can be found in the notebook in `notebooks/Data Processing.ipynb`.

## Generating and including Constellations
Satellite constellations (e.g. Starlink, OneWeb) are usually described by a list of orbital shells.
An orbital shell is described by a 4-tuple with information about `altitude`, `inclination`, `number of
planes`, and `number of satellites` per plane. We provide a notebook 
`notebooks/ConstellationGeneration/ConstellationGeneration.ipynb` that can be used 
to generate constellation data from orbital shell parameters.

### Generating constellation:
* Initialize the constellation by executing the first cell and providing metadata in the second cell (1)
* Create a shell by providing the 4 shell arguments, and further parameters (extra params) if necessary (2.1).
Store the temporary shell data by executing the cell (2.2)
* Turn satellites into position and velocity vectors by executing cell (3)
* Write the files by executing cell (4) and save them by executing cell (5)

A more detailed guide is included in the notebook.

### Including the constellation data in simulation (`io` section):

* In the configuration file for the simulation, include the constellation(s) by
defining `constellationList` and assigning the constellation name(s); Syntax: 
{name;}name ; e.g. Astra;Starlink;OneWeb
* Assign values for `constellationFrequency`, `constellationCutoff` and
`altitudeSpread` (more information below)

### How constellation satellites are inserted to the simulation

When starting the simulation with satellite constellations every
`constellationFrequency` iterations, satellites that are due to be added
are inserted to the simulation. They are only inserted if there is no other
object within the range `constellationCutoff` of the position of the new 
Satellite. If there are objects close enough to the new satellite, the insertion 
is postponed to the next time satellites are added. 

The insertion of a constellation takes as long as specified by the `duration` 
parameter in the respective .yaml file and happens shell by shell. The time it 
takes to insert one shell of a constellation depends on the percentage of 
satellites the shell contributes to. Satellites of each orbital 
shell are inserted plane by plane and linearly over time.

Satellites are inserted with varying altitude based on a random variable that is
normally distributed. `AltitudeSpread` determines within what range most satellites 
(99.74%) deviate from `altitude`. The parameter equals 3 times the standard deviation
of the normal distribution.





## Output

LADDS has multiple options for output that can be (de)activated mostly independent of each other via YAML. See `cfg/default_cfg.yaml` for relevant options.

### VTK
`.vtu` files in XML/ASCII layout that can be loaded into [Paraview](https://www.paraview.org/) for visualization.

### HDF5
A single `.h5` containing particle and conjunction data from a full simulation run with the following structure:
```
/
├── CollisionData
│   └── <IterationNr>
│       └── (Dataset) Collisions
│           idA idB distanceSquared
└── ParticleData
    └── <IterationNr>
        └── Particles
            ├── (Dataset) IDs
            ├── (Dataset) Positions
            │   x y z
            └── (Dataset) Velocities
                x y z
```

Collision data is tracked every iteration, particle data only in intervals that are defined in the YAML file. To keep file size reasonable compression is supported.

### CSV
If HDF5 output is disabled entirely, collision data is written in a `.csv` file in ASCII layout.


