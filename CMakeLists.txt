cmake_minimum_required(VERSION 2.6.0)

# pull in the pods macros. See cmake/pods.cmake for documentation
set(POD_NAME p2l-point-process-core)
include(cmake/pods.cmake)

# automatically build LCM types.  This also defines a number of CMake
# variables, see cmake/lcmtypes.cmake for details
#include(cmake/lcmtypes.cmake)
#lcmtypes_build()

# include_directories(${LCMTYPES_INCLUDE_DIRS})

#set(CMAKE_CXX_FLAGS "-std=c++0x -pedantic -Wall -O2 -g")
#set(CMAKE_CXX_FLAGS "-std=c++0x -pedantic -Wall -O0 -g3")

add_definitions( -std=c++0x -Wall -fdiagnostics-show-option )
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3")

option ( USE_PEDANTIC "Turn on -pendantic mode in gcc. This will spit out *lots* of warnings from lcm :-(, but hopefully none from the rest of the code" OFF)
if( USE_PEDANTIC )
    add_definitions( -pedantic )
endif (USE_PEDANTIC)


# the point-process-core library
add_library( p2l-point-process-core SHARED
  src/point_math.cpp
  src/entropy.cpp
  src/marked_grid.cpp
  src/context.cpp
  )
pods_install_headers( 
  src/point_math.hpp
  src/marked_grid.hpp
  src/entropy.hpp
  src/point_process.hpp
  src/context.hpp
  DESTINATION
  point-process-core )
pods_use_pkg_config_packages(p2l-point-process-core 
    gsl-1.16 
    boost-1.54.0 
    p2l-math-core
    lcmtypes_p2l-math-core
    p2l-probability-core
    lcmtypes_p2l-probability-core)
pods_install_libraries( p2l-point-process-core )
pods_install_pkg_config_file(p2l-point-process-core
    CFLAGS
    LIBS -lp2l-point-process-core
    REQUIRES gsl-1.16 boost-1.54.0 p2l-math-core lcmtypes_p2l-math-core p2l-probability-core lcmtypes_p2l-probability-core 
    VERSION 0.0.1)


add_subdirectory( test )
