project(NM CXX)
cmake_minimum_required(VERSION 3.6)

set(CMAKE_VERBOSE_MAKEFILE off)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -std=c++14")
set(NM_INCLUDES)
set(NM_LIBRARIES)

set (CMAKE_CXX_STANDARD 14)
cmake_policy(SET CMP0048 NEW)

# include
include_directories("${CMAKE_SOURCE_DIR}")
include_directories(external/plog/include/)

# GSL (external)
find_package(GSL)
if(GSL_FOUND)
    set(HAVE_LIBGSL 1)
    list(APPEND NM_INCLUDES ${GSL_INCLUDE_DIR})
    list(APPEND NM_LIBRARIES ${GSL_LIBRARIES})
else(GSL_FOUND)
    message(FATAL_ERROR "Require GSL: set GSL_HOME")
endif(GSL_FOUND)

# create the output dir
#file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# add source files
set(NM_SRCS
    src/gsl.cpp
    src/timer.cpp
    src/utils.cpp
)

# make library
include_directories (include ${NM_INCLUDES})
add_library( NM_LIB ${NM_SRCS} ${NM_INCLUDES} )

# make executables
# add_executable(trapezium_v0 src/integration/trapezium_v0.cpp)
# target_link_libraries (trapezium_v0 NM_LIB ${NM_LIBRARIES})

# add_executable(trapezium_v1 src/integration/trapezium_v1.cpp)
# target_link_libraries (trapezium_v1 NM_LIB ${NM_LIBRARIES})

# add_executable(trapezium_v2 src/integration/trapezium_v2.cpp)
# target_link_libraries (trapezium_v2 NM_LIB ${NM_LIBRARIES})

# add_executable(trapezium_v3 src/integration/trapezium_v3.cpp)
# target_link_libraries (trapezium_v3 NM_LIB ${NM_LIBRARIES})

# add_executable(montecarlo src/integration/montecarlo.cpp)
# target_link_libraries (montecarlo NM_LIB ${NM_LIBRARIES})

# add_executable(tridiag src/tridiag/tridiag.cpp)
# target_link_libraries (tridiag NM_LIB ${NM_LIBRARIES})

# add_executable(rk2 src/ode/rk2.cpp)
# target_link_libraries (rk2 NM_LIB ${NM_LIBRARIES})

# add_executable(lorenz src/ode/lorenz.cpp)
# target_link_libraries (lorenz NM_LIB ${NM_LIBRARIES})

# add_executable(advection src/pde/advection.cpp)
# target_link_libraries (advection NM_LIB ${NM_LIBRARIES})

# add_executable(advectionCN src/pde/advection_CN.cpp)
# target_link_libraries (advectionCN NM_LIB ${NM_LIBRARIES})

# add_executable(diffusion src/pde/diffusion.cpp)
# target_link_libraries (diffusion NM_LIB ${NM_LIBRARIES})

# add_executable(meanvalue src/mcmc/meanvalue.cpp)
# target_link_libraries (meanvalue NM_LIB ${NM_LIBRARIES})

# add_executable(isampling src/mcmc/isampling.cpp)
# target_link_libraries (isampling NM_LIB ${NM_LIBRARIES})

# add_executable(diagnosing src/misc/diagnosing.cpp)
# target_link_libraries (diagnosing NM_LIB ${NM_LIBRARIES})

add_executable(adiabatic src/examples/adiabatic.cpp)
target_link_libraries (adiabatic NM_LIB ${NM_LIBRARIES})

add_executable(galaxy src/examples/galaxy.cpp)
target_link_libraries (galaxy NM_LIB ${NM_LIBRARIES})

# show summary
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "GRYPHON_INCLUDES ...... = ${NM_INCLUDES}")
message (STATUS "GRYPHON_LIBRARIES ..... = ${NM_LIBRARIES}")
