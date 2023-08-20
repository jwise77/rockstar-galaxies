# From: https://git.ligo.org/cds/software/advligorts/
# FindRPC
# ---------
#
# Find the RPC includes and library.
#
#  - They can either be in the Transport-Independent RPC (TIRPC) package or 
#  still in glibc.    
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project:
#  ``RPC_FOUND``
#   true if the RPC header are found from TIRPC or glibc
#  ``RPC_INCLUDE_DIRS```
#   where to find rpc.h, etc. (could be simply /usr/include)
#  ``RPC_LIBRARIES``
#   the libraries to link against for RPC (may be empty if glibc)
# ``TIRPC_FOUND``
#   true if the TIRPC headers and libraries were found.
#
# The following may be set if TIRPC found
# ``TIRPC_INCLUDE_DIRS``
#   where to find rpc.h, etc. from TIRPC
# ``TIRPC_LIBRARIES``
#   the libraries to link against to use TIRPC.
# ``TIRPC_VERSION``
#   the version of TIRPC found.
#
if (RPC_FOUND)
else (RPC_FOUND)

  if (TIRPC_FOUND)
  else (TIRPC_FOUND)

    find_package(PkgConfig QUIET)
    pkg_check_modules(PC_TIRPC libtirpc)

    find_path(TIRPC_INCLUDE_DIRS
      NAMES netconfig.h
      PATH_SUFFIXES tirpc
      HINTS ${PC_TIRPC_INCLUDE_DIRS}
    )

    find_library(TIRPC_LIBRARIES
      NAMES tirpc
      HINTS ${PC_TIRPC_LIBRARY_DIRS}
    )

    set(TIRPC_VERSION ${PC_TIRPC_VERSION})

    include(FindPackageHandleStandardArgs)

    find_package_handle_standard_args(TIRPC
       REQUIRED_VARS TIRPC_LIBRARIES TIRPC_INCLUDE_DIRS
       VERSION_VAR TIRPC_VERSION
    )

    mark_as_advanced(TIRPC_INCLUDE_DIRS TIRPC_LIBRARIES)
 
  endif (TIRPC_FOUND)
  
  if (TIRPC_FOUND)
    set (RPC_FOUND TRUE)
    set (RPC_INCLUDE_DIRS ${TIRPC_INCLUDE_DIRS})
    set (RPC_LIBRARIES ${TIRPC_LIBRARIES})
  else (TIRPC_FOUND)
    find_path(rpc_RPC_H rpc.h hints /usr/include/rpc)
    if (rpc_RPC_H)
      Message("RPC found from glibc")
      set (RPC_FOUND TRUE)
      set (RPC_INCLUDE_DIRS "/usr/include")
      set (RPC_LIBRARIES "")
    else (rpc_RPC_H)
      Message("RPC Libraries not Found")
    endif (rpc_RPC_H)
  endif (TIRPC_FOUND)

endif (RPC_FOUND)
