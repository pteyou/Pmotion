cmake_minimum_required(VERSION 3.0)
project(pedestrianMotion)

add_subdirectory(FMM)

set(SourcesList4 "forces.cpp" "reunion.cpp")
add_executable(reunion ${SourcesList4})
target_link_libraries(reunion FMM)

set(SourcesList5 "forces.cpp" "cinema.cpp")
add_executable(cinema ${SourcesList5})
target_link_libraries(cinema FMM)

set(SourcesList6 "forces.cpp" "placePeopleInZenith.cpp" "zenith.cpp")
add_executable(zenith ${SourcesList6})
target_link_libraries(zenith FMM)

include_directories("~/Eigen3/")
include_directories("~/boost_1_69_0/")
set (CMAKE_CXX_FLAGS "-O3 -g")