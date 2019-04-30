cmake_minimum_required(VERSION 3.0)
project(FMM)


set(SourcesList1 "FMM.cpp" "rootFinding.cpp" "../utilities.cpp")
add_library(FMM ${SourcesList1})

include_directories("~/Eigen3/")
include_directories("~/boost_1_69_0/")
set (CMAKE_CXX_FLAGS "-O3 -g")