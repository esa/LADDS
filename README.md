# LADDS - Large-scale Deterministic Debris Simulation

Codebase for the ARIADNA Study between TU Munich and ESA's Advanced Concepts Team.

## Requirements
* CMake >= 3.14
* make (build-essentials, or equivalent)
* A C++17 Compiler (recommended: gcc >=7 or clang >=8)
* OpenMP >= 4.5
* git (for fetching dependencies)

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
TODO
