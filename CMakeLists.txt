cmake_minimum_required (VERSION 3.1)

SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++11 -pthread ")

project (COR-Predictor)

add_executable(COR-Predictor COR-main.cpp)
