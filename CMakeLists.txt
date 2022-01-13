cmake_minimum_required(VERSION 3.12)

# Top-level options
option(ENABLE_CUDA "Enable CUDA, which might be necessary for some of the plugins." OFF)

set(languages C CXX)
if (ENABLE_CUDA)
  list(APPEND languages CUDA)
endif()

# Create top-level project
project(VISUS_ParaView_Plugins ${languages})

# Set common include path
set(common_include ${CMAKE_CURRENT_SOURCE_DIR})

# Copy markdown readme files into documentation folder
file(GLOB_RECURSE md_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/plugins/*.md)

foreach(md_file ${md_files})
  get_filename_component(md_file_path ${md_file} DIRECTORY)

  string(REPLACE "plugins/" "" md_file_path "${md_file_path}")
  string(REPLACE "modules/" "" md_file_name "${md_file_path}")

  install(FILES ${md_file} DESTINATION share/paraview-plugins/doc RENAME ${md_file_name}.md)
endforeach()

install(FILES "README.md" DESTINATION share/paraview-plugins/doc RENAME Readme.md)

# Add plugins
add_subdirectory(plugins)
