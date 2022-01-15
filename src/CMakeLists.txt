cmake_minimum_required(VERSION 2.8)

project(hhta)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "-g -ggdb -Wall -Wno-unknown-pragmas -fPIC -DIL_STD -fexceptions -std=c++0x -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native -openmp -O3 -DNDEBUG")

set(Boost_NO_SYSTEM_PATHS OFF)
find_package(Boost COMPONENTS program_options timer system chrono filesystem)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

include_directories("${PROJECT_BINARY_DIR}")
add_executable(hhta hhta.cpp)
target_link_libraries(hhta ${Boost_LIBRARIES})