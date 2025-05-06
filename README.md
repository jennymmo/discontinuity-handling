# ODE-integrators-for-Lagrangian-particles

This repository contains the code used to run simulations and analyse results for a paper submitted to GMD. The simulation code is implemented in Fortran, using three external libraries (netCDF, hdf5 and bspline-fortran). The analysis and plotting of results is done in jupyter notebooks.

## Structure of the repo

### build

This directory is where you build your simulation code. It initially contains only a .gitignore file ignoring everything but itself, ensuring that files in this folder do not clutter the git status.

### data

This directory contain the ocean current datasets used to run the simulations, as well as the initial positions used for the 10000 particles. For reference, it also contains a map showing the extent of the datasets. The ocean current data are obtained from MET Norway, for additional details see data/README.md

### notebooks

This directory contains jupyter notebooks that are used to analyse and plot the results of running the simulations.

### results

This directory is where the simulation results will be stared. It initially contains only a .gitignore file ignoring everything but itself, ensuring that files in this folder do not clutter the git status.

### src

This directory contains the fortran source code to run the simulations.

## Build instructions for fortran code

* Download and build the bspline-fortran library from https://github.com/jacobwilliams/bspline-fortran
* Make sure you have
  * cmake version 3.5 or higher
  * hdf5, with fortran interface
  * netcdf, with fortran interface
  * a recent version of gfortran
* Go into the build directory (or create your own, somewhere)
* We use cmake to configure the build, and we need to pass information about the location of the different libraries. The following works for me on xubuntu 20.04 with netcdf installed from apt-get. Adjust the paths to fit your system.
* `cmake .. -DBSPLINE_LIBRARY_PATH=bspline_path -DNETCF_LIBRARY_PATH=netcdf_libs -DNETCDF_INCLUDE_DIRS=netcdf_include -DCMAKE_BUILD_TYPE=RELEASE`
  * `bspline_path` should be the folder where `bspline_module.o` and `bspline_module.mod` etc. are found. This depends on where you built bspline-fortran.
  * `netcdf_libs` should be the folder where `libnetcdf.so` and `libnetcdff.so` found. On my system, running xubuntu 18.04 with `libnetcdf-dev` and `libnetcdff-dev` installed from apt-get, these files are found in `/usr/lib/x86_64-linux-gnu`.
  * `netcdf_include` should be the folder where `netcdf.h` and `netcdf.mod` etc. are found. On my system, running xubuntu 18.04 with `libnetcdf-dev` and `libnetcdff-dev` installed from apt-get, these files are found in `/usr/include`.
* Then running `make` should build the project. If that doesn't work, try running `make VERBOSE=1`, copy the last command, where the error occured, and see if all the paths are correct.
* Running `make` will build four executables. How do use them is described below.

## Running the code

After building, you will have four executables, two for standard forwards runs, and two for the backtracking case. The forwards executable are `run_standard`, which runs a numerical experiment with the regular fixed-step integrators and `run_discontinuity_handling`, which runs an experiment with improved handling of the discontinuities. They both take the same arguments, indicating which dataset to use as input, and which order (order = polynomial degeer + 1) of interpolation to use. For example to run with the data from NorKyst 800, with linear interpolation, use

`  > ./run_normal norkyst800 2`

Run the program without any arguments to see the options.

A numerical experiment consists of running the code for:
* One dataset
* One order of interpolation
* Five different integrators
* 11 different timesteps

Some of the timesteps, for the reference solutions, are very short, and the simulations will take a long time to run. To change this, edit `run_standard.f90` and `run_discontinuity_handling.f90` as needed, the list of timesteps to run is found in the variable `timesteps` in both files.

For the backtracking case, run `run_backwards_standard` or `run_backwards_discontiunity_handling`. Both of these are set up to run the case presented in the paper, which uses the 4 km dataset, and linear interpolation.

## Instructions for using the jupyter notebooks

* First, build the fortran code and run the simulations.
* Then, stepping through the cells in order should reproduce the key results from the paper.
