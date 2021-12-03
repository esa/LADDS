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

## Processing TLE Input 
Data on current satellites etc. is often found [online](https://www.space-track.org/) in the [TLE format](https://en.wikipedia.org/wiki/Two-line_element_set). We include a Jupyter notebook which can be used to process TLE data with pykep to create and analyze suitable datasets. Detailed instructions can be found in the notebook in `notebooks/Data Processing.ipynb`.
