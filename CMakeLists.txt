cmake_minimum_required(VERSION 3.10)

#Set the project name
project(KERR_SIM VERSION 1.0 LANGUAGES CXX Fortran)

#Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#Add the executable
add_executable(KERR_SIM src/main.cpp src/ynogk.f90)

# Add option to use the accretion disk or not
option(USE_DISK "Use the accretion disk" ON)

#find packages
find_package(OpenGL REQUIRED)
find_package(glfw3 3.3 REQUIRED)
if(NOT glfw3_FOUND)
    message(FATAL_ERROR "GLFW not found!")
endif()
find_package(GLEW REQUIRED)
find_package(GLU REQUIRED)
find_package(GLUT REQUIRED)

#Add include directories
target_include_directories(KERR_SIM PUBLIC include)

#Include GLFW and GLEW Headers
target_include_directories(KERR_SIM PRIVATE ${GLFW_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS} ${OPENGL__INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} ${GLU_INCLUDE_DIRS})

#Define the paths to the shader files
target_compile_definitions(KERR_SIM PRIVATE
    VERTEX_SHADER_PATH="${CMAKE_SOURCE_DIR}/include/Kerr_Sim/Shaders/vShader.vs"
    FRAGMENT_SHADER_PATH="${CMAKE_SOURCE_DIR}/include/Kerr_Sim/Shaders/fShader.fs"
)

#Link GLFW and GLEW libraries
target_link_libraries(KERR_SIM PRIVATE glfw ${GLEW_LIBRARIES} ${OPENGL_LIBRARIES} ${GLU_LIBRARIES} ${GLUT_LIBRARIES} m)

#Debugging Information
message(STATUS "OpenGL libraries: ${OPENGL_LIBRARIES}")
message(STATUS "GLFW3_libraries: ${GLFW3_LIBRARIES}")
message(STATUS "GLEW libraries : ${GLEW_LIBRARIES}")
message(STATUS "GLUT libraries: ${freeglut3_libraries}")
message(STATUS "library directory: ${CMAKE_CURRENT_SOURCE_DIR}/libs")