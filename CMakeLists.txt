cmake_minimum_required(VERSION 3.20)

# project name and version
project(matrix VERSION 0.1.0)

configure_file(matrixConfig.h.in matrixConfig.h)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# add the matrix library
include_directories(includes)

add_subdirectory(test)


