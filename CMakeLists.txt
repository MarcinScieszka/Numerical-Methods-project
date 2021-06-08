cmake_minimum_required(VERSION 3.19)
project(mo_lab11)

set(CMAKE_CXX_STANDARD 20)

add_library(
        my_libraries
        calerf.h
        calerf.cpp
)

add_executable(mo_lab11 main.cpp)

target_link_libraries(mo_lab11 PRIVATE my_libraries)