add_executable(${PROJECT_NAME} dune-geometric-pde.cc)
target_link_dune_default_libraries(${PROJECT_NAME})

add_dune_alberta_flags(WORLDDIM 2 ${PROJECT_NAME})
add_dune_suitesparse_flags(${PROJECT_NAME})

add_definitions(-DSOURCEDIR="${PROJECT_SOURCE_DIR}/src")
add_definitions(-DMSHFILESDIR="${PROJECT_SOURCE_DIR}/msh-files")
add_definitions(-DGRIDDIM=ALBERTA_DIM-1)
add_definitions(-DWORLDDIM=ALBERTA_DIM)
