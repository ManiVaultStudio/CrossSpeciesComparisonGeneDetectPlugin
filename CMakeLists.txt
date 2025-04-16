cmake_minimum_required(VERSION 3.22)

option(MV_UNITY_BUILD "Combine target source files into batches for faster compilation" OFF)

# ----------------------------------------------------------------------------- 
# Project
# -----------------------------------------------------------------------------
project(SimianViewerPairExpCompPlugin)

# ----------------------------------------------------------------------------- 
# CMake Options
# -----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DWIN32 /EHsc /MP /permissive- /Zc:__cplusplus")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
endif()

# ----------------------------------------------------------------------------- 
# Dependencies
# -----------------------------------------------------------------------------
find_package(Qt6 REQUIRED COMPONENTS Widgets WebEngineWidgets Concurrent)
find_package(ManiVault COMPONENTS Core PointData ClusterData CONFIG QUIET)

if(NOT DEFINED MV_CSCTD_INSTALL_DIR)
    set(MV_CSCTD_INSTALL_DIR "${MV_INSTALL_DIR}" CACHE PATH "Directory where TreeData is installed")
endif()

# ----------------------------------------------------------------------------- 
# Source Files
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

set(AUX
    res/CrossSpeciesComparisonGeneDetect_chart.qrc
)

qt6_add_resources(RESOURCE_FILES ${AUX})

source_group("Plugin" FILES ${PLUGIN_SOURCES})
source_group("LibFiles" FILES ${LIBS})
source_group("ResourceFiles" FILES ${RESOURCE_FILES})
source_group("ResourceFiles\\generated" FILES "${CMAKE_CURRENT_BINARY_DIR}/qrc_gene_detect_resources.cpp")
source_group("Aux" FILES ${AUX})

# ----------------------------------------------------------------------------- 
# Target Definition
# -----------------------------------------------------------------------------
add_library(${PROJECT_NAME} SHARED 
    ${PLUGIN_SOURCES}
    ${LIBS}
    ${RESOURCE_FILES}
)

# ----------------------------------------------------------------------------- 
# Include Directories
# -----------------------------------------------------------------------------
target_include_directories(${PROJECT_NAME} PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Common"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Data"
    "${ManiVault_INCLUDE_DIR}"
    "${MV_CSCTD_INSTALL_DIR}/$<CONFIGURATION>/include/"
)

# ----------------------------------------------------------------------------- 
# Compiler Features & Properties
# -----------------------------------------------------------------------------
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)

set_target_properties(${PROJECT_NAME} PROPERTIES
    FOLDER "ViewPlugins"
)

if(MV_UNITY_BUILD)
    set_target_properties(${PROJECT_NAME} PROPERTIES UNITY_BUILD ON)
endif()

# ----------------------------------------------------------------------------- 
# Linking
# -----------------------------------------------------------------------------
target_link_libraries(${PROJECT_NAME} PRIVATE
    Qt6::Widgets
    Qt6::WebEngineWidgets
    Qt6::Concurrent
    ManiVault::Core
    ManiVault::PointData
    ManiVault::ClusterData
)

set(CROSSSPECIESCOMPARISONTREEDATA_LINK_LIBRARY
    "${CSCTDPLUGIN_LINK_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}CrossSpeciesComparisonTreeData${MV_LINK_SUFFIX}"
)

target_link_libraries(${PROJECT_NAME} PRIVATE "${CROSSSPECIESCOMPARISONTREEDATA_LINK_LIBRARY}")

# ----------------------------------------------------------------------------- 
# Installation
# -----------------------------------------------------------------------------
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION Plugins COMPONENT PLUGINS
    LIBRARY DESTINATION Plugins COMPONENT PLUGINS
)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND "${CMAKE_COMMAND}"
        --install "${CMAKE_CURRENT_BINARY_DIR}"
        --config "$<CONFIGURATION>"
        --prefix "${ManiVault_INSTALL_DIR}/$<CONFIGURATION>"
)

# ----------------------------------------------------------------------------- 
# MSVC Debug Settings
# -----------------------------------------------------------------------------
if(MSVC)
    set_property(TARGET ${PROJECT_NAME} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY 
        "$<IF:$<CONFIG:DEBUG>,${ManiVault_INSTALL_DIR}/debug,${ManiVault_INSTALL_DIR}/release>"
    )
    set_property(TARGET ${PROJECT_NAME} PROPERTY VS_DEBUGGER_COMMAND 
        "$<IF:$<CONFIG:DEBUG>,${ManiVault_INSTALL_DIR}/debug/ManiVault Studio.exe,${ManiVault_INSTALL_DIR}/release/ManiVault Studio.exe>"
    )
endif()
