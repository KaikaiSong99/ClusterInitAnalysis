cmake_minimum_required(VERSION 3.17)

# -----------------------------------------------------------------------------
# ClusterInitAnalysis Plugin
# -----------------------------------------------------------------------------
set(CLUSTERINITANALYSIS "ClusterInitAnalysisPlugin")
PROJECT(${CLUSTERINITANALYSIS})

# -----------------------------------------------------------------------------
# CMake Options
# -----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DWIN32 /EHsc /MP /permissive- /Zc:__cplusplus")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
endif(MSVC)

# -----------------------------------------------------------------------------
# Set install directory
# -----------------------------------------------------------------------------
# Check if the directory to the ManiVault installation has been provided
if(NOT DEFINED MV_INSTALL_DIR)
    set(MV_INSTALL_DIR "" CACHE PATH "Directory where ManiVault is installed")
    message(FATAL_ERROR "Please set MV_INSTALL_DIR to the directory where ManiVault is installed")
endif()
file(TO_CMAKE_PATH ${MV_INSTALL_DIR} MV_INSTALL_DIR)

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------
find_package(Qt6 COMPONENTS Widgets WebEngineWidgets REQUIRED)

# -----------------------------------------------------------------------------
# Source files
# -----------------------------------------------------------------------------
# Define the plugin sources
set(CLUSTERINITANALYSIS_SOURCES
    src/ClusterInitAnalysisPlugin.h
    src/ClusterInitAnalysisPlugin.cpp
    src/SettingsAction.h
    src/SettingsAction.cpp
    src/ClusterInitAnalysisPlugin.json
)

set(PLUGIN_MOC_HEADERS
	src/ClusterInitAnalysisPlugin.h
)

source_group( Plugin FILES ${CLUSTERINITANALYSIS_SOURCES})

# -----------------------------------------------------------------------------
# CMake Target
# -----------------------------------------------------------------------------
# Create dynamic library for the plugin
add_library(${CLUSTERINITANALYSIS} SHARED ${CLUSTERINITANALYSIS_SOURCES})

qt_wrap_cpp(CLUSTERINITANALYSIS_MOC ${PLUGIN_MOC_HEADERS} TARGET ${CLUSTERINITANALYSIS})
target_sources(${CLUSTERINITANALYSIS} PRIVATE ${CLUSTERINITANALYSIS_MOC})

# -----------------------------------------------------------------------------
# Target include directories
# -----------------------------------------------------------------------------
# Include ManiVault headers, including system data plugins
target_include_directories(${CLUSTERINITANALYSIS} PRIVATE "${MV_INSTALL_DIR}/$<CONFIGURATION>/include/")

# -----------------------------------------------------------------------------
# Target properties
# -----------------------------------------------------------------------------
# Request C++17
target_compile_features(${CLUSTERINITANALYSIS} PRIVATE cxx_std_17)

# Link to Qt libraries
target_link_libraries(${CLUSTERINITANALYSIS} PRIVATE Qt6::Widgets)
target_link_libraries(${CLUSTERINITANALYSIS} PRIVATE Qt6::WebEngineWidgets)

# -----------------------------------------------------------------------------
# Target library linking
# -----------------------------------------------------------------------------
# Link to ManiVault and data plugins
# The link path in this repo assumes that the ManiVault core was built locally
# in contrast to having been installed with an installer. In the latter case you'll have
# to adapt the MV_LINK_PATH and PLUGIN_LINK_PATH to your install folder
set(MV_LINK_PATH "${MV_INSTALL_DIR}/$<CONFIGURATION>/lib")
set(PLUGIN_LINK_PATH "${MV_INSTALL_DIR}/$<CONFIGURATION>/$<IF:$<CXX_COMPILER_ID:MSVC>,lib,Plugins>")
set(MV_LINK_SUFFIX $<IF:$<CXX_COMPILER_ID:MSVC>,${CMAKE_LINK_LIBRARY_SUFFIX},${CMAKE_SHARED_LIBRARY_SUFFIX}>)

set(MV_LINK_LIBRARY "${MV_LINK_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}MV_Public${MV_LINK_SUFFIX}")
set(POINTDATA_LINK_LIBRARY "${PLUGIN_LINK_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}PointData${MV_LINK_SUFFIX}") 
set(CLUSTERDATA_LINK_LIBRARY "${PLUGIN_LINK_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}ClusterData${MV_LINK_SUFFIX}") 

target_link_libraries(${CLUSTERINITANALYSIS} PRIVATE "${MV_LINK_LIBRARY}")
target_link_libraries(${CLUSTERINITANALYSIS} PRIVATE "${POINTDATA_LINK_LIBRARY}")
target_link_libraries(${CLUSTERINITANALYSIS} PRIVATE "${CLUSTERDATA_LINK_LIBRARY}")

# -----------------------------------------------------------------------------
# Target installation
# -----------------------------------------------------------------------------
# Install the shared plugin libary to the "Plugins" folder in the ManiVault install directory
install(TARGETS ${CLUSTERINITANALYSIS}
    RUNTIME DESTINATION Plugins COMPONENT PLUGINS # Windows .dll
    LIBRARY DESTINATION Plugins COMPONENT PLUGINS # Linux/Mac .so
)

add_custom_command(TARGET ${CLUSTERINITANALYSIS} POST_BUILD
    COMMAND "${CMAKE_COMMAND}"
        --install ${CMAKE_CURRENT_BINARY_DIR}
        --config $<CONFIGURATION>
        --prefix ${MV_INSTALL_DIR}/$<CONFIGURATION>
)

# -----------------------------------------------------------------------------
# Miscellaneous
# -----------------------------------------------------------------------------
# Automatically set the debug environment (command + working directory) for MSVC
if(MSVC)
    set_property(TARGET ${CLUSTERINITANALYSIS} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY $<IF:$<CONFIG:DEBUG>,${MV_INSTALL_DIR}/debug,${MV_INSTALL_DIR}/release>)
    set_property(TARGET ${CLUSTERINITANALYSIS} PROPERTY VS_DEBUGGER_COMMAND $<IF:$<CONFIG:DEBUG>,"${MV_INSTALL_DIR}/debug/ManiVault Studio.exe","${MV_INSTALL_DIR}/release/ManiVault Studio.exe">)
endif()
