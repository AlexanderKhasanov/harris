cmake_minimum_required(VERSION 3.8)
project(Harris)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp detector_harrisa.hpp detector_harrisa.cpp)

find_package( OpenCV REQUIRED )

add_executable(Harris ${SOURCE_FILES})

target_link_libraries( Harris ${OpenCV_LIBS} )