# Find MPI
find_package(MPI REQUIRED)
mark_as_advanced(MPI_EXTRA_LIBRARY)
mark_as_advanced(MPI_LIBRARY)
option (ENABLE_MPI "Enable the compilation of the MPI communication code" on)

