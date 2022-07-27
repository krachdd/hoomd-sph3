# the list of all source files
set(HIPPER_TEST_FILES
    test_hipper.cc
    test_kernel.cu
    test_cub.cu
)
# the list of test labels to turn into executables
set(HIPPER_TEST_LABELS
    kernel
    CUB
)

# look for catch quietly, and pull down v2.12.2 if not found
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()
find_package(Catch2 QUIET)
if(NOT Catch2_FOUND)
    message(STATUS "Downloading Catch2...")
    set(url https://github.com/catchorg/Catch2/releases/download/v2.12.2/catch.hpp)
    file(DOWNLOAD ${url} "${CMAKE_CURRENT_BINARY_DIR}/catch2/catch.hpp" STATUS status)
    list(GET status 0 error)
    if(error)
        message(FATAL_ERROR "... could not download from ${url}.")
    endif()
    set(Catch2_ROOT ${CMAKE_CURRENT_BINARY_DIR})
    find_package(Catch2 REQUIRED)
endif()

# set full path to files for use in subdirectories
list(TRANSFORM HIPPER_TEST_FILES PREPEND "${CMAKE_CURRENT_LIST_DIR}/")

add_subdirectory(cuda)
add_subdirectory(hip)
