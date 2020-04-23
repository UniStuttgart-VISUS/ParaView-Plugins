cmake_minimum_required(VERSION 3.12)

pv_module(attachment_separation_lines ${PROJECT_NAME} "" attachment_separation_lines_target)

# Find and link Eigen
find_package(Eigen3 REQUIRED NO_MODULE)

if(NOT TARGET Eigen3::Eigen)
  message(FATAL_ERROR "Target for eigen3 not found")
endif()

target_link_libraries(${attachment_separation_lines_target} PRIVATE Eigen3::Eigen)
