# Try to find the TIMEMORY library and headers
# Usage of this module is as follows
#
#     find_package( TIMEMORY )
#     if(TIMEMORY_FOUND)
#         include_directories(${TIMEMORY_INCLUDE_DIRS})
#         add_executable(foo foo.cc)
#         target_link_libraries(foo ${TIMEMORY_LIBRARIES})
#     endif()
#
# You can provide a minimum version number that should be used.
# If you provide this version number and specify the REQUIRED attribute,
# this module will fail if it can't find a TIMEMORY of the specified version
# or higher. If you further specify the EXACT attribute, then this module
# will fail if it can't find a TIMEMORY with a version eaxctly as specified.
#
# ===========================================================================
# Variables used by this module which can be used to change the default
# behaviour, and hence need to be set before calling find_package:
#
#  TIMEMORY_ROOT_DIR
#       The preferred installation prefix for searching for TIMEMORY
#       Set this if the module has problems finding the proper TIMEMORY installation.
#
# If you don't supply TIMEMORY_ROOT_DIR, the module will search on the standard
# system paths.
#
# ============================================================================
# Variables set by this module:
#
#  TIMEMORY_FOUND           System has TIMEMORY.
#
#  TIMEMORY_INCLUDE_DIRS    TIMEMORY include directories: not cached.
#
#  TIMEMORY_LIBRARIES       Link to these to use the TIMEMORY library: not cached.
#
# ===========================================================================
# If TIMEMORY is installed in a non-standard way, e.g. a non GNU-style install
# of <prefix>/{lib,include}, then this module may fail to locate the headers
# and libraries as needed. In this case, the following cached variables can
# be editted to point to the correct locations.
#
#  TIMEMORY_INCLUDE_DIR    The path to the TIMEMORY include directory: cached
#
#  TIMEMORY_LIBRARY        The path to the TIMEMORY library: cached
#
# You should not need to set these in the vast majority of cases
#

#----------------------------------------------------------------------------------------#

find_path(TIMEMORY_ROOT_DIR
    NAMES
        timemory/timemory.h timemory/timemory.hpp timemory/library.h
    HINTS
        ENV TIMEMORY_ROOT_DIR
        ENV CPATH
    PATH_SUFFIXES
        include bin lib lib64 
    DOC
        "TIMEMORY installation root directory")

#----------------------------------------------------------------------------------------#

find_path(TIMEMORY_INCLUDE_DIR
    NAMES
        timemory/timemory.h timemory/timemory.hpp timemory/library.h
    HINTS
        ${TIMEMORY_ROOT_DIR}
        ENV TIMEMORY_ROOT_DIR
        ENV CPATH
    PATH_SUFFIXES
        include/timemory timemory include 
    DOC
        "Path to the TIMEMORY headers")

#----------------------------------------------------------------------------------------#

find_library(TIMEMORY_LIBRARY
    NAMES
        timemory
    HINTS
        ${TIMEMORY_ROOT_DIR}
        ENV TIMEMORY_ROOT_DIR
        ENV LD_LIBRARY_PATH
        ENV LIBRARY_PATH
        ENV DYLD_LIBRARY_PATH
    PATH_SUFFIXES
        lib
        lib64
        timemory
        lib/timemory
        lib64/timemory
        # system processor
        lib/timemory/${CMAKE_SYSTEM_PROCESSOR}
        lib64/timemory/${CMAKE_SYSTEM_PROCESSOR}
        lib/timemory/${CMAKE_SYSTEM_PROCESSOR}/lib
        lib64/timemory/${CMAKE_SYSTEM_PROCESSOR}/lib64
        lib64/timemory/${CMAKE_SYSTEM_PROCESSOR}/lib
        ${CMAKE_SYSTEM_PROCESSOR}/lib
        ${CMAKE_SYSTEM_PROCESSOR}/lib/timemory
        ${CMAKE_SYSTEM_PROCESSOR}/lib64
        ${CMAKE_SYSTEM_PROCESSOR}/lib64/timemory
    DOC
        "Path to the TIMEMORY library")

#----------------------------------------------------------------------------------------#

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set TIMEMORY_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(TIMEMORY DEFAULT_MSG TIMEMORY_INCLUDE_DIR TIMEMORY_LIBRARY)

#----------------------------------------------------------------------------------------#

if(TIMEMORY_FOUND)
    add_library(TIMEMORY INTERFACE)
    target_link_libraries(TIMEMORY INTERFACE ${TIMEMORY_LIBRARY})
    target_include_directories(TIMEMORY INTERFACE ${TIMEMORY_INCLUDE_DIR})
    get_filename_component(TIMEMORY_INCLUDE_DIRS ${TIMEMORY_INCLUDE_DIR} REALPATH)
    get_filename_component(TIMEMORY_LIBRARIES ${TIMEMORY_LIBRARY} REALPATH)
endif()

mark_as_advanced(TIMEMORY_INCLUDE_DIR TIMEMORY_LIBRARY)
