DUFS
===

[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)]

`DUFS` is a C++ implementation of the Directed Unbiased Frontier Sampling method, which generalizes Frontier Sampling and the Directed Unbiased Random Walk.

Features
------

  - Example script run_DUFS.sh

Requirements
------

To install this project, please ensure that you have installed the following (install guides are provided on the respective websites):

  - [Git](http://git-scm.com)
  - A C++ compiler, e.g., [GCC](https://gcc.gnu.org/), [clang](http://clang.llvm.org/), [MinGW](http://www.mingw.org/)
  - [CMake](http://www.cmake.org "CMake homepage")
  - [GSL](https://www.gnu.org/software/gsl/)


Installation
------


Project structure
-------------

This project has been set up with a specific file/folder structure in mind. The following describes some important features of this setup:

  - `cmake/Modules` : Contains `CMake` modules
  - `include`: Project header files (*.hpp)
  - `examples`: Shell scripts
  - `src`: Project source files (*.cpp)
  - `CMakeLists.txt`: main `CMakelists.txt` file for project
  - `Dependencies.cmake`: list of dependencies and automated build, triggered if dependency cannot be found locally
  - `LICENSE.md`: license file for project
  - `ProjectFiles.cmake`: list of project source files to build
  - `README.md`: project readme file


TODO
------------
