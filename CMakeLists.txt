cmake_minimum_required(VERSION 3.8)
project(Harris)

find_package(Boost 1.57.0 COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp detector_harrisa.hpp detector_harrisa.cpp)

find_package( OpenCV REQUIRED )

add_executable(Harris ${SOURCE_FILES})

target_link_libraries( Harris ${OpenCV_LIBS} )
target_link_libraries( Harris ${Boost_LIBRARIES} )