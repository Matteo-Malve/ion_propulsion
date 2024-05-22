cmake_minimum_required(VERSION 3.22)

find_package(deal.II 9.3.3 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
deal_ii_initialize_cached_variables()

project(ion-propulsion)

#--------------------------------------------------

include_directories(src)

#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/results)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/results)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)

set(SOURCE_FILES
    src/main.cpp
)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
deal_ii_setup_target(${PROJECT_NAME})

#--------------------------------------------------


add_custom_target(run
  COMMAND ${PROJECT_NAME}
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/../results
  )