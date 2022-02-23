# LADDS - Large-scale Deterministic Debris Simulation

Codebase for the ARIADNA Study between TU Munich and ESA's Advanced Concepts Team. A more detailed project description can be found on the [Advanced Concept Team's webpage](https://www.esa.int/gsp/ACT/projects/debris_hpc/).

## Requirements
* CMake >= 3.19
* make (build-essentials, or equivalent)
* A C++17 Compiler (recommended: gcc >=7 or clang >=8, only clang-10 is tested)
* OpenMP >= 4.5
* git (for fetching dependencies)
* [TBB](https://github.com/oneapi-src/oneTBB) (Breakup-Model needs this)

### Optional
* Doxygen
* clang-format-9

## Important Dependencies
The following codes play an important role in this project. They are downloaded and managed via CMake at configure time:
* [AutoPas](https://github.com/AutoPas/AutoPas)
* [NASA Breakup Model](https://github.com/esa/NASA-breakup-model-cpp)
* [Orbit Propagator](https://github.com/FG-TUM/OrbitPropagator)

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

## Running
The simulation requires one `yaml` file as argument which specifies the necessary options.
```bash
./ladds myInput.yaml
```
For an overview of all possible options see [`cfg/default_cfg.yaml`](cfg/default_cfg.yaml). Most parameters have
a default value which is used when they are left unspecified. The full configuration, including defaulted
values, is shown in the console output when executing the simulation.

### Checkpoints
LADDS features a checkpoint mechanic where a simulation can be restarted from an HDF5 file.
To create a checkpoint the simulation just needs to write an HDF5 file.
This is done via the following part of the configuration:
```yaml
io:
  csv:
    fileName: initial_population.csv  # = input
  hdf5:
    fileName: checkpointFile.h5       # = output
    writeFrequency: 100               # how often LADDS writes to HDF5  
```

The next simulation, which starts from `checkpointFile.h5` then needs to have the following:
```yaml
io:
  hdf5:
    checkpoint:
      file: checkpointFile.h5         # = input AND output
      iteration: 500                  # Iteration where to start from.
    writeFrequency: 100               # how often LADDS writes to HDF5  
```
LADDS will append any new data to the checkpoint file.

**IMPORTANT**: All file paths are relative to the [`data`](data/) directory!

## Simulating Breakups
The code is capable to simulate fatal collisions between two bodies via the 
[NASA Breakup Model](https://github.com/esa/NASA-breakup-model-cpp). This feature can be activated in
the `yaml` file via:
```yaml
sim:
  breakup:
    enabled: true
```

## Calibrating AutoPas
By default, some robust but static configuration is set by LADDS. You can change this by specifying the
algorithmic options AutoPas is allowed to use in the YAML file. If more than one configuration can be built
from these options AutoPas will tune over them at run time.

### Enable Auto Tuning
If you are unsure what algorithmic configuration you want to use for AutoPas just let AutoPas guide you.
For this, the following needs be activated in the `yaml` file:
```yaml
autopas:
  tuningMode: true
```
In this mode, the simulation is only executed for one AutoPas tuning-phase. At the end of this phase, a copy
of the full configuration is created, which contains the algorithm configuration that AutoPas deemed to be
the fastest. This configuration can then be used to run the actual simulation at optimal speed.

### Analyzing AutoPas Configurations
AutoPas can be compiled to dump information about the performance of the algorithms it uses to `.csv` files. 
For this set the CMake variables:  `AUTOPAS_LOG_TUNINGDATA=ON` and `AUTOPAS_LOG_TUNINGRESULTS=ON`.
* `AutoPas_tuningData.csv` contains the timing data of all samples AutoPas collected.
* `AutoPas_tuningResults.csv` contains the result of each tuning phase.

## Processing TLE Input 
Data on current satellites etc. is often found [online](https://www.space-track.org/) in the [TLE format](https://en.wikipedia.org/wiki/Two-line_element_set). We include a Jupyter notebook which can be used to process TLE data with pykep to create and analyze suitable datasets. Detailed instructions can be found in the notebook in `notebooks/Data Processing.ipynb`.

## Output

LADDS has multiple options for output that can be (de)activated mostly independent of each other via YAML. See [`cfg/default_cfg.yaml`](cfg/default_cfg.yaml) for relevant options.

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
    │   └── Particles
    │       ├── (Dataset) IDs
    │       ├── (Dataset) Positions
    │       │   x y z
    │       └── (Dataset) Velocities
    │           x y z
    └─── (Dataset) ConstantProperties
         id cosparId mass radius bcInv activityState
```

Collision data is tracked every iteration, particle data only in intervals that are defined in the YAML file.
`ConstantProperties` contains properties of all particles that existed over the course of the simulation. 
Due to burn ups or breakups the number of particles in any iteration might differ but `id`s are unique! 
To keep file size reasonable compression is supported.

### CSV
If HDF5 output is disabled entirely, collision data is written in a `.csv` file in ASCII layout.

