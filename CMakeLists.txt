cmake_minimum_required(VERSION 3.6)
project(Mètodes_Numèrics_2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    basic_structs.h
    ex14.c)
set(EX15_FILES
        basic_structs.h
        ex15.c)
set(EX21_FILES
        ex21.c)

add_executable(ex14 ${SOURCE_FILES})
add_executable(ex15 ${EX15_FILES})
add_executable(ex21 ${EX21_FILES})

target_link_libraries(ex21 m)