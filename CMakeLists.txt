cmake_minimum_required(VERSION 3.14)
project(unscentedKalmanFIlter)

set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(include)

add_library(ukf src/unscented_kalman_filter.cpp)
add_executable(unscentedKalmanFilter main.cpp)
target_link_libraries(unscentedKalmanFilter ukf)
