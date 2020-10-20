cmake_minimum_required(VERSION 3.10)
project(PtDAlgorithms)

set(CMAKE_CXX_STANDARD 98)

set(CMAKE_BINARY_DIR .)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS_DEBUG "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
add_executable(cpp_test src/test.cpp src/phase.h src/generation.h src/phase.cpp src/generation.cpp src/io.cpp src/io.h)
add_library(PtDAlgorithms src/phase.cpp src/phase.h src/io.cpp src/io.h)
add_executable(tester2 src/tester.cpp src/phase.cpp src/phase.h src/io.cpp src/io.h)

