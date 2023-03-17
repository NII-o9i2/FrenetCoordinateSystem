cmake_minimum_required(VERSION 3.24)
project(frenet_coordinate_project)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-O2)

#add_subdirectory(osqp-cpp)
find_package(osqp REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(benchmark REQUIRED)
#include_directories(utils/qp_solver)

add_executable(frenet_coordinate_project
        frenet_coordinate_system/frenet_coordinate_main.cpp
        utils/math_utils.cc
        utils/spline.cpp
        utils/vec2d.cc
        utils/line_segment2d.cc
        utils/cartesian_frenet_conversion.cc
        frenet_coordinate_system/frenet_coordinate_system_new.cpp
        frenet_coordinate_system/frenet_coordinate_system.cpp
        main.cpp
        utils/qp_solver/osqp_solver.cpp)
target_include_directories(frenet_coordinate_project PUBLIC
        fplus/
        utils/
)
target_link_libraries(frenet_coordinate_project benchmark osqp)



