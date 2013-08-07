cmake_minimum_required(VERSION 2.6.0)

# pull in the pods macros. See cmake/pods.cmake for documentation
set(POD_NAME point-process-core)
include(cmake/pods.cmake)

# automatically build LCM types.  This also defines a number of CMake
# variables, see cmake/lcmtypes.cmake for details
include(cmake/lcmtypes.cmake)
lcmtypes_build()

include_directories(${LCMTYPES_INCLUDE_DIRS})

#set(CMAKE_CXX_FLAGS "-std=c++0x -pedantic -Wall -O2 -g")
#set(CMAKE_CXX_FLAGS "-std=c++0x -pedantic -Wall -O0 -g3")

add_definitions( -std=c++0x -pedantic -Wall )
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3")


# the point-process-core library
add_library( point-process-core SHARED
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
pods_use_pkg_config_packages(point-process-core gsl boost-1.51.0 lcmtypes_point-process-core)
pods_install_libraries( point-process-core )
pods_install_pkg_config_file(point-process-core
    CFLAGS
    LIBS -lpoint-process-core
    REQUIRES gsl boost-1.51.0
    VERSION 0.0.1)


add_subdirectory( test )
