cmake_minimum_required(VERSION 3.5)
project (VariableTimestep)
enable_language (Fortran)


# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set (CMAKE_Fortran_FLAGS "-Wall")
set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -pg -fimplicit-none -fbacktrace -fcheck=all -fbounds-check -ffpe-trap=zero")
set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -mtune=native")


# Files with modules and stuff
set (fortran_modules
    "currentdata.f90"
    "experiment.f90"
    "input.f90"
    "integrator.f90"
    "interpolator.f90"
    "output.f90"
    "parameters.f90"
    )


# Using hdf5 to store output
#find_package (HDF5 COMPONENTS HL Fortran_HL Fortran)

# Using netCDF to read current data from files
file (GLOB  NETCDF_LIBRARIES  "${NETCDF_LIBRARY_PATH}/libnetcdf*.so")
file (GLOB  HDF5_LIBRARIES  "${NETCDF_LIBRARY_PATH}/libhdf*.so")
include_directories (${HDF5_INCLUDE_DIRS})
include_directories (${NETCDF_INCLUDE_DIRS})
message ("++ HDF5_LIBRARIES:       ${HDF5_LIBRARIES}")
message ("++ HDF5_INCLUDE_DIRS:    ${HDF5_INCLUDE_DIRS}")
message ("++ NETCDF_LIBRARIES:     ${NETCDF_LIBRARIES}")
message ("++ NETCDF_INCLUDE_DIRS:  ${NETCDF_INCLUDE_DIRS}")

# Need to link to the bspline-fortran library
# https://github.com/jacobwilliams/bspline-fortran
file (GLOB BSPLINE_LIBRARIES "${BSPLINE_LIBRARY_PATH}/build/lib/*.a")
include_directories (${BSPLINE_LIBRARY_PATH}/build)
include_directories (${BSPLINE_LIBRARY_PATH}/build/lib)
message ("++ BSPLINE_LIBRARY_PATH: ${BSPLINE_LIBRARY_PATH}")
message ("++ BSPLINE_LIBRARIES:    ${BSPLINE_LIBRARIES}")

#### Program to run all timesteps for one dataset and one order
add_executable (run_normal "run_normal.f90" ${fortran_modules})
target_link_libraries (run_normal ${BSPLINE_LIBRARIES})
target_link_libraries (run_normal ${HDF5_LIBRARIES})
target_link_libraries (run_normal ${NETCDF_LIBRARIES})

#### Program to run all timesteps for one dataset and one order
add_executable (run_discontinuity_handling "run_discontinuity_handling.f90" ${fortran_modules})
target_link_libraries (run_discontinuity_handling ${BSPLINE_LIBRARIES})
target_link_libraries (run_discontinuity_handling ${HDF5_LIBRARIES})
target_link_libraries (run_discontinuity_handling ${NETCDF_LIBRARIES})
