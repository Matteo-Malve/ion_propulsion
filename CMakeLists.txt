cmake_minimum_required(VERSION 3.22)

#--------------------------------------------------------------------
# Set default build type if not specified (Debug or Release)
# This only works for single-configuration generators (e.g. Makefiles)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build (Debug or Release)" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
#--------------------------------------------------------------------

# Find deal.II (adjust the path as needed)
find_package(deal.II 9.3.3 REQUIRED
        /Applications/deal.II.app/Contents/Resources/spack/opt/dealii-9.6.0-rc1-3efk/lib/cmake/deal.II
)
deal_ii_initialize_cached_variables()

# Define your target name
set(TARGET "ion_propulsion")

# Set the project name (using a dash is fine, but note that the executable name will match ${PROJECT_NAME})
project(ion-propulsion)

#--------------------------------------------------------------------
# Include directories and output settings
include_directories(src include)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../results)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../results)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Specify your source files
set(SOURCE_FILES
        src/Evaluation.cpp
        src/Framework.cpp
        src/Data.cpp
        src/ion_propulsion.cpp
        src/DualFunctional.cpp
        src/LaplaceSolver.cpp
        src/Refinement.cpp
        src/Constants.cpp
)

# (Optional) OS X SDK specification
set(CMAKE_OSX_SYSROOT "/Library/Developer/CommandLineTools/SDKs/MacOSX14.4.sdk")

# This variable isn’t used by the manual target setup, but is left for reference:
set(TARGET_SRC src/${TARGET}.cpp)
#--------------------------------------------------------------------

# Create the executable from all source files
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
deal_ii_setup_target(${PROJECT_NAME})

#--------------------------------------------------------------------
# Custom target to run the program after building it
add_custom_target(run
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --parallel
        COMMAND ${PROJECT_NAME}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/../results
)

#--------------------------------------------------------------------
# Custom target for "make distclean"
# This target removes CMake-generated build files and the results directory.
# Note: Running "distclean" means you will need to re-run CMake to reconfigure your build.
add_custom_target(distclean
        COMMAND ${CMAKE_COMMAND} -E echo "Performing distclean..."
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/CMakeFiles
        COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/CMakeCache.txt
        COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/cmake_install.cmake
        COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/compile_commands.json
        COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/Makefile
        COMMENT "Removing all generated build files and results directory."
)
