cmake_minimum_required(VERSION 3.13.4)

find_package(deal.II 9.3.3 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
deal_ii_initialize_cached_variables()

project(second-grid)


#--------------------------------------------------

add_executable(second_grid second_grid.cpp)
deal_ii_setup_target(second_grid)

#--------------------------------------------------


add_custom_target(run
  COMMAND second_grid
  )
