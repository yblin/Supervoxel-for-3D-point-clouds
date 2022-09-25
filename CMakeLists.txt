cmake_minimum_required(VERSION 3.20)
project(supervoxel)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

file(GLOB_RECURSE HEADER_FILES
    codelibrary
)

file(GLOB_RECURSE SOURCE_FILES
    codelibrary/*.cpp
)

add_executable(supervoxel main.cc ${SOURCE_FILES})
target_include_directories(supervoxel PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")