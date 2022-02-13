cmake_minimum_required(VERSION 3.10)
project(PtDAlgorithms)

set(CMAKE_CXX_STANDARD 98)

set(CMAKE_BINARY_DIR .)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
#set(CMAKE_CXX_FLAGS_DEBUG "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

cmake_minimum_required(VERSION 3.9)
project(ptdalgorithms VERSION 1.0.0 DESCRIPTION "Efficient graph based phase-type distribution algorithms")
include(GNUInstallDirs)
add_library(libptdalgorithms SHARED api/c/ptdalgorithms.h src/c/ptdalgorithms.c)
set_target_properties(libptdalgorithms PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION 1
        PUBLIC_HEADER api/c/ptdalgorithms.h)
configure_file(ptdalgorithms.pc.in ptdalgorithms.pc @ONLY)
target_include_directories(libptdalgorithms PRIVATE src/c)
install(TARGETS libptdalgorithms
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_BINARY_DIR}/ptdalgorithms.pc
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

add_library(libptdalgorithmscpp SHARED api/c/ptdalgorithms.h src/c/ptdalgorithms.c src/c/ptdalgorithms.h src/cpp/ptdalgorithmscpp.cpp api/cpp/ptdalgorithmscpp.h src/cpp/ptdalgorithmscpp.h)

set_target_properties(libptdalgorithmscpp PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION 1
        PUBLIC_HEADER api/cpp/ptdalgorithmscpp.h)
configure_file(ptdalgorithmscpp.pc.in ptdalgorithmscpp.pc @ONLY)
target_include_directories(libptdalgorithmscpp PRIVATE src/cpp)
install(TARGETS libptdalgorithmscpp
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_BINARY_DIR}/ptdalgorithmscpp.pc
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)


#add_executable(testcpp test/testcpp.cpp)
#target_link_libraries (testcpp LINK_PUBLIC libptdalgorithmscpp)

#add_executable(tester test/tester.cpp)
#target_link_libraries (tester LINK_PUBLIC libptdalgorithms)

