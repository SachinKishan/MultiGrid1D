# Minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Project name
project(MG VERSION 1.0 LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Source files
set(SOURCES
    source.cpp
    imgui/imgui.cpp
    imgui/imgui_demo.cpp
    imgui/imgui_draw.cpp
    imgui/imgui_tables.cpp
    imgui/imgui_widgets.cpp
    imgui/imgui_impl_glfw.cpp
    imgui/imgui_impl_opengl3.cpp
    imgui/implot.cpp
    imgui/implot_demo.cpp
    imgui/implot_items.cpp
)

# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/include/GLFW
    ${PROJECT_SOURCE_DIR}/include/glm
    ${PROJECT_SOURCE_DIR}/include/KHR
    ${PROJECT_SOURCE_DIR}/imgui
    ${PROJECT_SOURCE_DIR}/Eigen
)

# Add the GLAD library
add_subdirectory(vendor/glad)

# Add executable
add_executable(MG ${SOURCES})

# Platform-specific configurations
if(WIN32)
    # Windows-specific configuration
    target_link_libraries(MG PRIVATE glad ${PROJECT_SOURCE_DIR}/libs/glfw3.lib)
elseif(APPLE)
    # macOS-specific configuration
    find_package(PkgConfig REQUIRED)
    pkg_search_module(GLFW REQUIRED glfw3)
    include_directories(${GLFW_INCLUDE_DIRS})
    target_link_libraries(MG PRIVATE glad ${GLFW_LIBRARIES} "-framework OpenGL" "-framework Cocoa" "-framework IOKit" "-framework CoreVideo")
endif()



# Build types
set(CMAKE_BUILD_TYPE Debug)

