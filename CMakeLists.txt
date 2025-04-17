cmake_minimum_required(VERSION 3.22)

option(MV_UNITY_BUILD "Combine target source files into batches for faster compilation" OFF)

# -----------------------------------------------------------------------------
# Project Setup
# -----------------------------------------------------------------------------
project(CrossSpeciesComparisonGeneDetectPlugin LANGUAGES CXX)

# -----------------------------------------------------------------------------
# CMake Options
# -----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

if(MSVC)
    add_compile_options(
        /DWIN32
        /EHsc
        /MP
        /permissive-
        /Zc:__cplusplus
    )
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
endif()

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------
find_package(Qt6 COMPONENTS Widgets WebEngineWidgets Concurrent REQUIRED)
find_package(ManiVault COMPONENTS Core PointData ClusterData CONFIG REQUIRED)

# Find CrossSpeciesComparisonTreeData plugin
if(NOT DEFINED MV_CSCTD_INSTALL_DIR)
    set(MV_CSCTD_INSTALL_DIR ${MV_INSTALL_DIR} CACHE PATH "Directory where TreeData is installed")
    message(STATUS "MV_CSCTD_INSTALL_DIR not defined. Using MV_INSTALL_DIR as default: ${MV_CSCTD_INSTALL_DIR}")
endif()

# Set paths for TreeData plugin
set(CSCTD_INCLUDE_DIR "${MV_CSCTD_INSTALL_DIR}/$<CONFIG>/include/CrossSpeciesComparisonTreeData")
set(CSCTD_LIB_DIR "${MV_CSCTD_INSTALL_DIR}/$<CONFIG>/Plugins")
set(CSCTD_LIB_NAME "CrossSpeciesComparisonTreeData")

# -----------------------------------------------------------------------------
# Source files
# -----------------------------------------------------------------------------
set(PLUGIN_SOURCES
    src/CrossSpeciesComparisonGeneDetectPlugin.h
    src/CrossSpeciesComparisonGeneDetectPlugin.cpp
    src/SettingsAction.h
    src/SettingsAction.cpp
    src/CrossSpeciesComparisonGeneDetectPlugin.json
)

set(LIBS
    src/lib/Clustering/fastcluster.cpp
    src/lib/Clustering/fastcluster.h
    src/lib/Clustering/fastcluster_dm.h
    src/lib/Clustering/fastcluster_R_dm.h
    src/lib/JSONnlohmann/json.hpp
    src/lib/JSONnlohmann/json_fwd.hpp
    src/lib/Distance/annoylib.h
    src/lib/Distance/kissrandom.h
    src/lib/Distance/mman.h
    src/lib/NewickComparator/newick_comparator.h
)

set(PLUGIN_MOC_HEADERS
    src/CrossSpeciesComparisonGeneDetectPlugin.h
)

qt6_add_resources(RESOURCE_FILES res/CrossSpeciesComparisonGeneDetect_chart.qrc)

# Organize files in IDE
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${PLUGIN_SOURCES} ${LIBS} ${RESOURCE_FILES})

# -----------------------------------------------------------------------------
# Target Definition
# -----------------------------------------------------------------------------
add_library(${PROJECT_NAME} SHARED 
    ${PLUGIN_SOURCES} 
    ${LIBS} 
    ${RESOURCE_FILES}
)

# -----------------------------------------------------------------------------
# Target Properties
# -----------------------------------------------------------------------------
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "ViewPlugins")

if(MV_UNITY_BUILD)
    set_target_properties(${PROJECT_NAME} PROPERTIES UNITY_BUILD ON)
endif()

# -----------------------------------------------------------------------------
# Include Directories
# -----------------------------------------------------------------------------
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Common
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Data
    ${ManiVault_INCLUDE_DIR}
    ${CSCTD_INCLUDE_DIR}
)

# -----------------------------------------------------------------------------
# Link Libraries
# -----------------------------------------------------------------------------
target_link_libraries(${PROJECT_NAME} PRIVATE
    Qt6::Widgets
    Qt6::WebEngineWidgets
    Qt6::Concurrent
    ManiVault::Core
    ManiVault::PointData
    ManiVault::ClusterData
)

# Link against CrossSpeciesComparisonTreeData
target_link_directories(${PROJECT_NAME} PRIVATE ${CSCTD_LIB_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE ${CSCTD_LIB_NAME})

# -----------------------------------------------------------------------------
# Installation
# -----------------------------------------------------------------------------
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION Plugins COMPONENT PLUGINS
    LIBRARY DESTINATION Plugins COMPONENT PLUGINS
    ARCHIVE DESTINATION lib COMPONENT PLUGINS
)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND}
        --install ${CMAKE_CURRENT_BINARY_DIR}
        --config $<CONFIG>
        --prefix ${ManiVault_INSTALL_DIR}/$<CONFIG>
)

# -----------------------------------------------------------------------------
# MSVC Debug Configuration
# -----------------------------------------------------------------------------
if(MSVC)
    set_property(TARGET ${PROJECT_NAME} PROPERTY 
        VS_DEBUGGER_WORKING_DIRECTORY "${ManiVault_INSTALL_DIR}/$<CONFIG>")
    set_property(TARGET ${PROJECT_NAME} PROPERTY 
        VS_DEBUGGER_COMMAND "${ManiVault_INSTALL_DIR}/$<CONFIG>/ManiVault Studio.exe")
endif()
