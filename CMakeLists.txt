cmake_minimum_required(VERSION 2.8.3)
project(singularity_4points)

set(CMAKE_BUILD_TYPE Release)

find_package(VISP REQUIRED COMPONENTS core vs imgproc mbt robot)

set(BUILD_EXAMPLES OFF)
set(ENABLE_DYNAMIC_PLOT  OFF)
set(BUILD_PARSER ON)
add_subdirectory(submodules/log2plot)

set(BUILD_TESTS OFF)   
add_subdirectory(submodules/opengv)
add_subdirectory(submodules/lambdatwist_pnp)


add_definitions( -fPIC )
add_definitions (-march=native)
add_definitions(-DSRC_PATH="${CMAKE_CURRENT_SOURCE_DIR}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

include_directories(
  include
  ${VISP_INCLUDE_DIRS}
  submodules/log2plot/include
  submodules/opengv/include
  submodules/lambdatwist_pnp
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
target_link_libraries(${PROJECT_NAME} opengv lambdatwist ${VISP_LIBRARIES} yaml-cpp log2plot)

add_executable(rose  src/rose.cpp)
target_link_libraries(rose ${PROJECT_NAME})

add_executable(visual_servo  src/visual_servo.cpp)
target_link_libraries(visual_servo ${PROJECT_NAME})

add_executable(sphere src/sphere.cpp)
target_link_libraries(sphere  ${PROJECT_NAME})

add_executable(energy src/energy.cpp)
target_link_libraries(energy ${PROJECT_NAME})

add_executable(best_orientation_demo src/best_orientation_demo.cpp)
target_link_libraries(best_orientation_demo ${PROJECT_NAME})

