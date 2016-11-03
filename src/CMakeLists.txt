add_executable(${PROJECT_NAME} dune-geometric-pde.cc)
target_link_dune_default_libraries(${PROJECT_NAME})

add_definitions(-DSOURCEDIR="${PROJECT_SOURCE_DIR}")
add_definitions(-DMSHFILESDIR="$ENV{HOME}/mshfiles")

add_dune_alberta_flags(WORLDDIM 2 ${PROJECT_NAME})
add_dune_suitesparse_flags(${PROJECT_NAME})
