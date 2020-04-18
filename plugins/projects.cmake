cmake_minimum_required(VERSION 3.12)

function(pv_project NAME)
  project(${NAME} CXX)

  if (ParaView_VERSION VERSION_LESS 5.7)
    # Build modules
    add_subdirectory(modules)

    # Build plugin
    add_subdirectory(plugin)
  else()
    # Scan and build modules
    file(GLOB_RECURSE module_files modules/*.module)

    vtk_module_scan(
      MODULE_FILES                "${module_files}"
      PROVIDES_MODULES            modules
      WANT_BY_DEFAULT             ON
      HIDE_MODULES_FROM_CACHE     OFF
    )

    foreach(module ${modules})
      message(STATUS "  ...module '${module}'")
    endforeach()

    vtk_module_build(
      MODULES                     ${modules}
      INSTALL_HEADERS             OFF
    )

    # Scan and build plugins
    paraview_plugin_scan(
      PLUGIN_FILES                "${CMAKE_CURRENT_SOURCE_DIR}/plugin/${PROJECT_NAME}.plugin"
      PROVIDES_PLUGINS            plugins
      ENABLE_BY_DEFAULT           ON
      HIDE_PLUGINS_FROM_CACHE     OFF
    )

    message(STATUS "  ...plugin '${PROJECT_NAME}'")

    paraview_plugin_build(
      PLUGINS                     ${plugins}
    )
  endif()
endfunction()

function(pv_plugin NAME MODULES)
  set(module_names)
  foreach(module ${MODULES})
    list(APPEND module_names ${NAME}::${module})
  endforeach()

  if (ParaView_VERSION VERSION_LESS 5.7)
    message(STATUS "  ...plugin '${NAME}'")

    add_paraview_plugin(${NAME} 1.0
      SERVER_MANAGER_XML      ${NAME}.xml
    )

    target_link_libraries(${NAME} PRIVATE ${module_names})

    install(TARGETS ${NAME})
  else()
    paraview_add_plugin(${NAME}
      VERSION                  1.0
      SERVER_MANAGER_XML      ${NAME}.xml
      MODULES                 ${module_names}
    )
  endif()
endfunction()

function(pv_module NAME PLUGIN SOURCES)
  if (ParaView_VERSION VERSION_LESS 5.7)
    add_library(${PLUGIN}_${NAME} SHARED ${NAME}.h ${NAME}.cxx ${SOURCES})
    add_library(${PLUGIN}::${NAME} ALIAS ${PLUGIN}_${NAME})

    target_link_libraries(${PLUGIN}_${NAME} PRIVATE ${VTK_LIBRARIES})
  else()
    vtk_module_add_module(${PLUGIN}::${NAME}
      SOURCES
        ${NAME}.cxx
      HEADERS
        ${NAME}.h
    )
  endif()
endfunction()
