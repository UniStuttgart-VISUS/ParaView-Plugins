cmake_minimum_required(VERSION 3.12)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Find ParaView
include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(BUILD_SHARED_LIBS ON)

find_package(ParaView REQUIRED)

if (ParaView_VERSION VERSION_LESS 5.7)
  include(${PARAVIEW_USE_FILE})
  include(ParaViewPlugins)
endif()

# Get project functions
include(projects.cmake)

# Add all subdirectories
file(GLOB children RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *)

foreach(child ${children})
  if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${child})
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${child}/CMakeLists.txt)
      message(STATUS "Adding plugin directory '${child}'...")
      add_subdirectory(${child})
    endif()
  endif()
endforeach()
