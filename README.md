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
│           idA, idB, distanceSquared
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

