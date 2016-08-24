cmake_minimum_required(VERSION 2.8.3)
project(pso)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
)

include_directories(
  ${catkin_INCLUDE_DIRS}
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

catkin_package()
add_executable(RF src/RF.cpp)
target_link_libraries(RF ${catkin_LIBRARIES})
add_dependencies(RF pso_gencpp)
