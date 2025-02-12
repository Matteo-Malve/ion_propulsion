cmake_minimum_required(VERSION 3.22)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build (Debug or Release)" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Find deal.II (adjust the path as needed)
find_package(deal.II 9.3.3 REQUIRED /usr/local)
deal_ii_initialize_cached_variables()

# Define your target name
set(TARGET "ion_propulsion")

# Set the project name
project(ion-propulsion)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

# Create the executable from all source files
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
deal_ii_setup_target(${PROJECT_NAME})

# Make sure the target uses C++17
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

# Custom target to run the program after building it
add_custom_target(run
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --parallel
        COMMAND ${PROJECT_NAME}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/../results
)

# Custom target for "make distclean"
add_custom_target(distclean
        COMMAND ${CMAKE_COMMAND} -E echo "Performing distclean..."
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/CMakeFiles
        COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/CMakeCache.txt
        COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/cmake_install.cmake
        COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/compile_commands.json
        COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/Makefile
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_NAME}
        COMMENT "Removing all generated build files and results directory."
)