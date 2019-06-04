cmake_minimum_required(VERSION 2.8.12)
project(pylibseq)

find_package(pybind11)
message(STATUS "Found pybind11: ${pybind11_VERSION}")
if(${pybind11_VERSION} VERSION_LESS '2.2.3')
    message(FATAL_ERROR "pybind11 version must be >= '2.2.3'")
endif()

option(USE_WEFFCPP "Use -Weffc++ during compilation" ON)
option(BUILD_UNIT_TESTS "Build C++ modules for unit tests" ON)

if (USE_WEFFCPP)
    add_compile_options(-Weffc++)
endif()

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")

include_directories(BEFORE libsequence/src/libsequence)
add_subdirectory(libsequence)

# if(BUILD_UNIT_TESTS)
#     add_subdirectory(tests)
# endif(BUILD_UNIT_TESTS)
