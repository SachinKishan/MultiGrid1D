# vendor/glad/CMakeLists.txt

cmake_minimum_required(VERSION 3.10)
project(glad)

# Add the GLAD library
add_library(glad STATIC src/glad.c)

# Include directories for GLAD
target_include_directories(glad PUBLIC include)