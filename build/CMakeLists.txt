cmake_minimum_required(VERSION 3.13.4)

find_package(deal.II 9.4.0 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
deal_ii_initialize_cached_variables()

project(ion-propulsion)

#--------------------------------------------------

add_subdirectory(${CMAKE_SOURCE_DIR}/../src ${CMAKE_SOURCE_DIR}/../src/build)

#--------------------------------------------------

add_custom_target(cleanres
  COMMAND rm *.svg *.vtu
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/../results
  )