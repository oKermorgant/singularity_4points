cmake_minimum_required(VERSION 2.8.3)
project(singularity_4points)

set(CMAKE_BUILD_TYPE Release)

find_package(log2plot)
find_package(VISP REQUIRED COMPONENTS core vs imgproc mbt robot)
find_package(opengv)
find_package(lambdatwist)

set(CMAKE_CXX_STANDARD 17)

add_definitions( -fPIC )
add_definitions (-march=native)

add_definitions(-DSRC_PATH="${CMAKE_CURRENT_SOURCE_DIR}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

## Specify additional locations of header files
## Your package locations should be listed before other locations
#include_directories(include)
# TODO: Check names of system library include directories (visp)
include_directories(
  include
  ${VISP_INCLUDE_DIRS}
  ${LOG2PLOT_INCLUDE_DIRS}
  ${opengv_INCLUDE_DIRS} ${LAMBDATWIST_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME} SHARED
  include/task.h
  include/scene.h
  src/scene.cpp
  include/pose_estim.h
  src/pose_estim.cpp
  include/statistics.h
  src/statistics.cpp
  )
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES} opengv ${LAMBDATWIST_LIBRARIES} ${VISP_LIBRARIES} yaml-cpp ${LOG2PLOT_LIBRARIES})

add_executable(rose  src/rose.cpp)
target_link_libraries(rose ${PROJECT_NAME})

add_executable(visual_servo  src/visual_servo.cpp)
target_link_libraries(visual_servo ${PROJECT_NAME})

add_executable(sphere src/sphere.cpp)
target_link_libraries(sphere  ${PROJECT_NAME})

add_executable(best_orientation_demo src/best_orientation_demo.cpp)
target_link_libraries(best_orientation_demo ${PROJECT_NAME})

add_executable(test_sphere test_sphere.cpp)
target_link_libraries(test_sphere ${PROJECT_NAME})
