![CMake](https://github.com/TUD-STKS/VocalTractLabBackend-dev/actions/workflows/cmake.yml/badge.svg) ![MSBuild](https://github.com/TUD-STKS/VocalTractLabBackend-dev/actions/workflows/msbuild.yml/badge.svg) [![Python examples](https://github.com/TUD-STKS/VocalTractLabBackend-dev/actions/workflows/python_examples.yml/badge.svg)](https://github.com/TUD-STKS/VocalTractLabBackend-dev/actions/workflows/python_examples.yml) [![MATLAB examples](https://github.com/TUD-STKS/VocalTractLabBackend-dev/actions/workflows/matlab_examples.yml/badge.svg)](https://github.com/TUD-STKS/VocalTractLabBackend-dev/actions/workflows/matlab_examples.yml) 

# VocalTractLabBackend-dev
This repo contains the VocalTractLab backend source code and the C/C++ API for on-going development work. *This is not the place for official stable releases!* 

You can find the bleeding edge builds here including new and experimental features, some of which may break your old projects. For official stable releases, please check [the VocalTractLab website](https://www.vocaltractlab.de).

Please feel free to fork this repo and make your own contributions, though! 

The main branch of this repo is reviewed on a semi-regular basis for inclusion into the official release.

This repo may be included in other repos as a submodule wherever the backend source code or the C/C++ API is needed.

## What you can find in this repo
This repo contains sources of the backend of the articulatory synthesizer [VocalTractLab](https://www.vocaltractlab.de). The backend can be accessed through a C/C++ API offering convenient C-style functions to provide the most commonly used functionality (e.g. converting a gestural score file into an audio file using default synthesis settings), or as a static C++ library for full access to all objects in the backend. This repo therefore contains separate projects/targets for both the C/C++ API and the C++ static library.

## Getting started
- Clone the current main branch:
```
git clone https://github.com/TUD-STKS/VocalTractLabBackend-dev
```

### Build using CMake (Windows, Linux, macOS)
- Get the latet release of [CMake for your platform](https://cmake.org/)
- Create a folder ``out`` inside the cloned repository folder
- Open a shell/command prompt and navigate to ``out``
- Configure the project and generate a build system:
```
cmake .. -DCMAKE_BUILD_TYPE=Release
```
- Build the library (still from within the folder ``out``)
```
cmake --build . --config Release
```
This will build both the API and the static library, as well as some unit testing executables. If you are only interested in one of those, you can specify the `--target` parameter:
```
cmake --build . --config Release --target VocalTractLabApi
```
or
```
cmake --build . --config Release --target VocalTractLabBackend
```
The final binary files are placed into the subdirectory ``lib`` (of the repository root). Do not look for them in the ``out`` directory, which is only used as a temporary folder for the build process and can be safely deleted once the binaries are built.

### Run tests using shell/command prompt
- Open a shell/command prompt
- Navigate to cloned repository folder (not the ``out`` folder from the build process)
- run the tests by executing ``VtlApiTests`` in the ``lib`` folder created by the build process:
```
./lib/VtlApiTests
```

### Build using Visual Studio 2019 (Windows)
- Open ``VocalTractLabApi.sln`` in the folder `build/msw`
- Build the project ``VocalTractLabApi`` or ``VocalTractLabBackend``

### Use the VocalTractLab backend in your own projects
To include the VocalTractLab backend into your own projects, add the folder `include` from this repository to your project's include directories and link against the API or static backend library in the folder `lib`. 

You can then include the API functions like so:
```cpp
#include "VocalTractLabApi/VocalTractLabApi.h"
```

C++-Objects can be included from the static backend library through their respective header:
```cpp
#include "VocalTractLabBackend/Speaker.h"  // or substitute your desired header here
```

## How to use the VocalTractLab backend
The VocalTractLab API functions are documented with extensive comments in `VocalTractLabApi.h`. Unfortunately, since they were not previously public, the classes beyond the API in the static library are not (yet) documented consistently.
