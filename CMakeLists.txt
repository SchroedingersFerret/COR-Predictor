cmake_minimum_required (VERSION 3.1)


if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -pthread ")
else ()
  set (CMAKE_CXX_STANDARD 11)
endif ()

project (COR-Predictor)

add_executable(COR-Predictor COR-main.cpp)
