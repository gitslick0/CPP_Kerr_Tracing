cmake_minimum_required(VERSION 3.10)

# Set the project name and version
project(YourProject VERSION 1.0 LANGUAGES CXX Fortran)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add src directory
# Optionally, create a build directory for Project1
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_SOURCE_DIR}/build/OpenGL")
# Invoke CMake for OpenGL Project
add_subdirectory(src/OpenGL ${CMAKE_SOURCE_DIR}/build/OpenGL)

# add Testing framework
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_SOURCE_DIR}/build/extern/googletest")
add_subdirectory(extern/googletest ${CMAKE_SOURCE_DIR}/build/extern/googletest)

#add_test(NAME OpenGL_TEST COMMAND ${CMAKE_SOURCE_DIR}/build/OpenGL/runTests)
enable_testing()