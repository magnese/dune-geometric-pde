find_package(Alberta REQUIRED)
include(AddAlbertaFlags)

find_package(SuiteSparse OPTIONAL_COMPONENTS LDL SPQR UMFPACK REQUIRED)
include(AddSuiteSparseFlags)
