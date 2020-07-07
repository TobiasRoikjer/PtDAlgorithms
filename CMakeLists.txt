cmake_minimum_required(VERSION 3.10)
project(PtDAlgorithms)

set(CMAKE_CXX_STANDARD 11)

add_executable(cpp_test src/test.cpp src/phase.cpp src/phase.h)
add_library(PtDAlgorithms src/phase.cpp src/phase.h)