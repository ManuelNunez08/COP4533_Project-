cmake_minimum_required(VERSION 3.21)
project(COP4533_Project)

# 'lib' is the folder with Google Test sources
add_subdirectory(lib)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_subdirectory(Testing)

set(CMAKE_CXX_STANDARD 14)

add_executable(COP4533_Project main.cpp)

add_executable(Google_Tests_run Testing/ProjectTests.cpp)
target_link_libraries(Google_Tests_run gtest gtest_main)