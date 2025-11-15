cmake_minimum_required(VERSION 3.10)

project(MAIN)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(include)
aux_source_directory(src DIR_SRCS)


add_executable(MAIN main.cpp ${DIR_SRCS})
