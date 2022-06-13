cmake_minimum_required(VERSION 3.12)

pv_module(swirling_helix ${PROJECT_NAME} "" swirling_helix_target)

# Find and link Eigen
find_package(Eigen3 REQUIRED NO_MODULE)

if(NOT TARGET Eigen3::Eigen)
  message(FATAL_ERROR "Target for eigen3 not found")
endif()

target_link_libraries(${swirling_helix_target} PUBLIC Eigen3::Eigen)
