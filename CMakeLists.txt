cmake_minimum_required(VERSION 3.10)

# Set the project name and version
project(KerrTracing VERSION 1.0 LANGUAGES CXX Fortran)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# CMake Module Directory
# Add the FindCfitsio
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/cmake/Modules/FindCFITSIO.cmake")
    message(WARNING "Didn't find CFITSIO.cmake file in Modules path")
endif()

# Set configuration
configure_file(configuration/root_directory.h.in configuration/root_directory.h)
include_directories(${CMAKE_BINARY_DIR}/configuration)

# Create relevant static libraries required for other projects
add_library(STB_IMAGE "src/stb_image.cpp")
target_include_directories(STB_IMAGE PUBLIC ${CMAKE_SOURCE_DIR}/include)
set(LIBS ${LIBS} STB_IMAGE)


# Optionally, create a build directory for Project1
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_SOURCE_DIR}/build/OpenGL")
# Invoke CMake for OpenGL Project
add_subdirectory(src/OpenGL ${CMAKE_SOURCE_DIR}/build/OpenGL)

# add Testing framework
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_SOURCE_DIR}/build/extern/googletest")
add_subdirectory(extern/googletest ${CMAKE_SOURCE_DIR}/build/extern/googletest)

# add Fitswrite chapter
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_SOURCE_DIR}/build/FITS_Write")
add_subdirectory(src/FITS_WRITE ${CMAKE_SOURCE_DIR}/build/FITS_Write)

#add_test(NAME OpenGL_TEST COMMAND ${CMAKE_SOURCE_DIR}/build/OpenGL/runTests)
enable_testing()