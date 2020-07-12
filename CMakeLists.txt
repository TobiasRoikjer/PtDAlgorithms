cmake_minimum_required(VERSION 3.10)
project(PtDAlgorithms)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_BINARY_DIR .)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(cpp_test src/test.cpp src/phase.cpp src/phase.h)
add_library(PtDAlgorithms src/phase.cpp src/phase.h)