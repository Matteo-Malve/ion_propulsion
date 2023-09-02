include_directories(Foundamentals GoalOrientedEstimator Mesh)


add_library(HelperFunctions OBJECT Foundamentals/HelperFunctions.cpp)
deal_ii_setup_target(HelperFunctions)

add_library(GridForge OBJECT Mesh/GridForge.cpp)
deal_ii_setup_target(GridForge)

add_library(Framework OBJECT Foundamentals/Framework.cpp)
deal_ii_setup_target(Framework)

add_library(ErrorController OBJECT Foundamentals/ErrorController.cpp)
deal_ii_setup_target(ErrorController)


add_library(Solver OBJECT Foundamentals/Solver.cpp)
deal_ii_setup_target(Solver)

add_library(PrimalSolver OBJECT Foundamentals/PrimalSolver.cpp)
deal_ii_setup_target(PrimalSolver)

add_library(DualFunctional OBJECT Foundamentals/DualFunctional.cpp)
deal_ii_setup_target(DualFunctional)

add_executable(main main.cpp $<TARGET_OBJECTS:HelperFunctions> $<TARGET_OBJECTS:GridForge> $<TARGET_OBJECTS:Framework> $<TARGET_OBJECTS:ErrorController> $<TARGET_OBJECTS:Solver> $<TARGET_OBJECTS:PrimalSolver> $<TARGET_OBJECTS:DualFunctional>)
deal_ii_setup_target(main)


#--------------------------------------------------

set_property(TARGET main
  PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../../results
  )
add_custom_target(run
  COMMAND main
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/../results
  )
