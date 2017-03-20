# Copyright Dmitry Anisimov danston@ymail.com (c) 2017.

# Open terminal
# cd path_to_the_folder/2d/bin/
# cmake -DCMAKE_BUILD_TYPE=Debug .. or cmake -DCMAKE_BUILD_TYPE=Release ..
# make
# ./blwp

cmake_minimum_required(VERSION 3.3)

project(blwp)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wextra -std=c++11 ")

set(SRC 
        extra/BarycentricCoordinatesR2.hpp
        extra/SegmentCoordinatesR2.hpp
        extra/VertexExpressionsR2.hpp
        extra/VertexR2.hpp
        coords/MeanValueR2.hpp
        coords/BlendedWachspressR2.hpp
        main.cpp
)

add_executable(blwp ${SRC})
