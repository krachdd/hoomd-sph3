cmake_minimum_required(VERSION 2.8.10...3.13 FATAL_ERROR)
project(pgsd)

set(PYTHON_MODULE_BASE_DIR "pgsd")

add_subdirectory(CMake)
include(CFlagsSetup)
include(LinterSetup)
include(MPISetup)


#############################################################################################
# Find libraries
include(PythonSetup)
include_directories(${PYTHON_INCLUDE_DIR})

if (WIN32)
add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

#############################################################################################
### Add project subdirs
include_directories(pgsd)
add_subdirectory(pgsd)
add_subdirectory(scripts)
