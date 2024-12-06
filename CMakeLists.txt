cmake_minimum_required(VERSION 3.22)

find_package(deal.II 9.3.3 REQUIRED
        /Applications/deal.II.app/Contents/Resources/spack/opt/dealii-9.6.0-rc1-3efk/lib/cmake/deal.II
)
deal_ii_initialize_cached_variables()

set(TARGET "ion_propulsion")

#--------------------------------------------------

include_directories(src include)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../results)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../results)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(SOURCE_FILES
    src/Evaluation.cpp
    src/Framework.cpp
    src/Data.cpp
    src/ion_propulsion.cpp
    src/DualFunctional.cpp
    src/LaplaceSolver.cpp
    src/Refinement.cpp
)
set(CMAKE_OSX_SYSROOT "/Library/Developer/CommandLineTools/SDKs/MacOSX14.4.sdk")

set(TARGET_SRC
    src/${TARGET}.cpp
)

project(ion-propulsion)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
deal_ii_setup_target(${PROJECT_NAME})

#--------------------------------------------------

# Custom target to run the program
add_custom_target(run
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --parallel
        COMMAND ${PROJECT_NAME}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/../results
)
