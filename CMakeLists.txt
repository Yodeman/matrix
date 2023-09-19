cmake_minimum_required(VERSION 3.20)

# project name and version
project(matrix VERSION 0.1.0)

configure_file(matrixConfig.h.in matrixConfig.h)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if (MSVC)
    add_compile_options(/W4 /EHsc)
else()
    add_compile_options(-Wall -Werror -Wpedantic)
endif()

# add the matrix library
include_directories(include)

# add_library(libmatrix OBJECT src/matrix.cpp)

add_subdirectory(test)


