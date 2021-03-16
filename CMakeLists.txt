project(NM CXX)
cmake_minimum_required(VERSION 3.6)

set(CMAKE_VERBOSE_MAKEFILE off)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-10")
set(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -std=c++11")
set(NM_INCLUDES)
set(NM_LIBRARIES)

set (CMAKE_CXX_STANDARD 11)
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
    src/timer.cpp
   )

# make library
include_directories (include ${NM_INCLUDES})
add_library( NM_LIB ${NM_SRCS} ${NM_INCLUDES} )

# make executables
add_executable(trapezium_v0 src/integration/trapezium_v0.cpp)
target_link_libraries (trapezium_v0 NM_LIB ${NM_LIBRARIES})

add_executable(trapezium_v1 src/integration/trapezium_v1.cpp)
target_link_libraries (trapezium_v1 NM_LIB ${NM_LIBRARIES})

add_executable(trapezium_v2 src/integration/trapezium_v2.cpp)
target_link_libraries (trapezium_v2 NM_LIB ${NM_LIBRARIES})

add_executable(trapezium_v3 src/integration/trapezium_v3.cpp)
target_link_libraries (trapezium_v3 NM_LIB ${NM_LIBRARIES})

add_executable(tridiag src/tridiag/tridiag.cpp)
target_link_libraries (tridiag NM_LIB ${NM_LIBRARIES})

# show summary
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "GRYPHON_INCLUDES ...... = ${NM_INCLUDES}")
message (STATUS "GRYPHON_LIBRARIES ..... = ${NM_LIBRARIES}")
