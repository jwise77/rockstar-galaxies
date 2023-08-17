CMAKE_MINIMUM_REQUIRED(VERSION 3.12)
project(
  ROCKSTAR
  LANGUAGES C
)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(/usr/include/tirpc)


if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif()

# Set C standard
set(CMAKE_C_STANDARD 99)

# Detect hdf5
find_package(HDF5 REQUIRED COMPONENTS C)
if (HDF5_FOUND)
  add_compile_definitions(H5_USE_16_API ENABLE_HDF5)
  link_libraries(HDF5::HDF5)
endif()


# add_compile_definitions(H5_USE_16_API ENABLE_HDF5)
add_compile_definitions(_LARGEFILE_SOURCE _LARGEFILE64_SOURCE _FILE_OFFSET_BITS=64 _DEFAULT_SOURCE _POSIX_SOURCE _POSIX_C_SOURCE=200809L _DARWIN_C_SOURCE)
add_compile_options(-fPIC -fno-math-errno -Wall -m64 -Og)
add_link_options(-fPIC)

# Create executable
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/io")
include_directories("${PROJECT_SOURCE_DIR}/inet")

project(ROCKSTAR C)

set(SOURCES
  main.c 
  rockstar.c 
  check_syscalls.c 
  fof.c 
  groupies.c 
  subhalo_metric.c 
  potential.c 
  nfw.c 
  jacobi.c 
  fun_times.c 
  interleaving.c 
  universe_time.c 
  hubble.c 
  integrate.c 
  distance.c 
  config_vars.c 
  config.c 
  bounds.c 
  inthash.c 
  client.c 
  server.c 
  merger.c 
  io/read_config.c 
  inet/socket.c 
  inet/rsocket.c 
  inet/address.c 
  io/meta_io.c 
  io/io_internal.c 
  io/io_ascii.c 
  io/stringparse.c 
  io/io_gadget.c 
  io/io_generic.c 
  io/io_art.c 
  io/io_nchilada.c 
  io/io_tipsy.c 
  io/io_bgc2.c 
  io/io_util.c 
  io/io_arepo.c 
  io/io_hdf5.c 
  io/io_enzo.c 
)

# Generate "rockstar-galaxies" executable
add_executable(rockstar-galaxies ${SOURCES})
target_link_libraries(rockstar-galaxies tirpc m)

# Generate "librockstar-galaxies.so" shared library
add_library(lib SHARED ${SOURCES})
target_link_libraries(lib tirpc m)
set_target_properties(lib PROPERTIES OUTPUT_NAME "rockstar-galaxies")

install(TARGETS rockstar-galaxies lib)
