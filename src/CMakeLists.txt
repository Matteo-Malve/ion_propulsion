include_directories(Foundamentals GoalOrientedEstimator Mesh)

add_library(Problem OBJECT Foundamentals/Problem.cpp)
deal_ii_setup_target(Problem)

add_library(HelperFunctions OBJECT Foundamentals/HelperFunctions.cpp)
deal_ii_setup_target(HelperFunctions)

add_library(GridForge OBJECT Mesh/GridForge.cpp)
deal_ii_setup_target(GridForge)

add_executable(main main.cpp $<TARGET_OBJECTS:Problem> $<TARGET_OBJECTS:HelperFunctions> $<TARGET_OBJECTS:GridForge>)
deal_ii_setup_target(main)

#--------------------------------------------------

set_property(TARGET main
  PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../../results
  )
add_custom_target(run
  COMMAND main
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/../results
  )
