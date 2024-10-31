cmake_minimum_required(VERSION 3.22)

find_package(deal.II 9.3.3 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
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
    src/problem.cpp
    src/main.cpp
)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
deal_ii_setup_target(${PROJECT_NAME})

#--------------------------------------------------

# Find Python library
find_package(PythonLibs 3 REQUIRED)

# Include Python directories
include_directories(${PYTHON_INCLUDE_DIRS})

# Link Python library
target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES})

# Custom target to run the program
add_custom_target(run
  COMMAND ${PROJECT_NAME}
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/../results
  )