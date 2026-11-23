cmake_minimum_required(VERSION 3.21)
project(test_package LANGUAGES CXX)

find_package(cnpy REQUIRED CONFIG)

add_executable(test_package test_package.cpp)
target_link_libraries(test_package PRIVATE cnpy::cnpy)
