# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
  set(FAR_CCD_TOPLEVEL_PROJECT OFF)
else()
  set(FAR_CCD_TOPLEVEL_PROJECT ON)
endif()

# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.18.0")
if(FAR_CCD_TOPLEVEL_PROJECT)
  cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
else()
  # Don't use cmake_minimum_required here to avoid implicitly overriding parent policies
  if(${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
    message(FATAL_ERROR "CMake required version to build FAR CCD is ${REQUIRED_CMAKE_VERSION}")
  endif()
endif()

# Include user-provided default options if available. We do that before the main
# `project()` so that we can define the C/C++ compilers from the option file.
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/FAR_CCD_Options.cmake)
  message(STATUS "Using local options file: ${CMAKE_CURRENT_SOURCE_DIR}/FAR_CCD_Options.cmake")
  include(${CMAKE_CURRENT_SOURCE_DIR}/FAR_CCD_Options.cmake)
endif()

# Enable ccache if available
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  option(FAR_CCD_WITH_CCACHE "Enable ccache when building FAR CCD" ${FAR_CCD_TOPLEVEL_PROJECT})
else()
  option(FAR_CCD_WITH_CCACHE "Enable ccache when building FAR CCD" OFF)
endif()
if(FAR_CCD_WITH_CCACHE AND CCACHE_PROGRAM)
  message(STATUS "Enabling Ccache support")
  set(ccacheEnv
    CCACHE_BASEDIR=${CMAKE_BINARY_DIR}
    CCACHE_SLOPPINESS=clang_index_store,include_file_ctime,include_file_mtime,locale,pch_defines,time_macros
  )
  foreach(lang IN ITEMS C CXX)
    set(CMAKE_${lang}_COMPILER_LAUNCHER
      ${CMAKE_COMMAND} -E env ${ccacheEnv} ${CCACHE_PROGRAM}
    )
  endforeach()
endif()


################################################################################
# CMake Policies
################################################################################

cmake_policy(SET CMP0054 NEW) # Only interpret if() arguments as variables or keywords when unquoted.
cmake_policy(SET CMP0076 NEW) # target_sources() command converts relative paths to absolute.
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24")
  cmake_policy(SET CMP0135 NEW) # Set the timestamps of all extracted contents to the time of the extraction.
endif()

################################################################################

project(FAR_CCD
        DESCRIPTION "Fast Root Approximate in Cubic CCD from [Lan et al. 2022]"
        LANGUAGES CXX
        VERSION "0.0.1")

option(FAR_CCD_BUILD_TESTS "Build unit-tests" ${FAR_CCD_TOPLEVEL_PROJECT})

# Set default minimum C++ standard
if(FAR_CCD_TOPLEVEL_PROJECT)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

### Configuration
set(FAR_CCD_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
set(FAR_CCD_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/far_ccd/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")

# General CMake utils
include(far_ccd_cpm_cache)
include(far_ccd_use_colors)

# Generate position-independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################
# FAR CCD Library
################################################################################

# Add an empty library and fill in the list of sources in `src/CMakeLists.txt`.
add_library(far_ccd)
add_library(far_ccd::far_ccd ALIAS far_ccd)

# Fill in configuration options
configure_file(
  "${FAR_CCD_SOURCE_DIR}/config.hpp.in"
  "${FAR_CCD_SOURCE_DIR}/config.hpp")

# Add source and header files to far_ccd
add_subdirectory("${FAR_CCD_SOURCE_DIR}")

# Public include directory for FAR CCD
target_include_directories(far_ccd PUBLIC "${FAR_CCD_INCLUDE_DIR}")

################################################################################
# Optional Definitions
################################################################################

# For MSVC, do not use the min and max macros.
target_compile_definitions(far_ccd PUBLIC NOMINMAX)

################################################################################
# Dependencies
################################################################################

# libigl
include(eigen)
target_link_libraries(far_ccd PUBLIC Eigen3::Eigen)

# Extra warnings (link last for highest priority)
include(far_ccd_warnings)
target_link_libraries(far_ccd PRIVATE far_ccd::warnings)

################################################################################
# Compiler options
################################################################################

# Use C++17
target_compile_features(far_ccd PUBLIC cxx_std_17)

################################################################################
# Tests
################################################################################

# Enable unit testing at the root level
if(FAR_CCD_TOPLEVEL_PROJECT AND FAR_CCD_BUILD_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()
