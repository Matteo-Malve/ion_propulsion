cmake_minimum_required(VERSION 3.22)

find_package(deal.II 9.3.3 REQUIRED
  /Applications/deal.II.app/Contents/Resources/spack/opt/dealii-9.6.0-rc1-3efk/lib/cmake/deal.II
  )
deal_ii_initialize_cached_variables()

project(ion-propulsion)

#--------------------------------------------------

include_directories(src, include)

#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/results)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/results)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)

set(SOURCE_FILES
    src/global.cpp
    src/utilities.cpp
    src/dual_functional.cpp
    src/problem_solver.cpp
    src/problem_error_refinement.cpp
    src/problem_convergence.cpp
    src/main.cpp
)
set(CMAKE_OSX_SYSROOT "/Library/Developer/CommandLineTools/SDKs/MacOSX14.4.sdk")


add_executable(${PROJECT_NAME} ${SOURCE_FILES})
deal_ii_setup_target(${PROJECT_NAME})

#--------------------------------------------------

# Find Python library
cmake_minimum_required(VERSION 3.12) # Ensure modern CMake
find_package(Python 3 REQUIRED COMPONENTS Interpreter Development)

include_directories(${Python_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${Python_LIBRARIES})

# Custom target to run the program
add_custom_target(run
  COMMAND ${PROJECT_NAME}
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/../results
  )