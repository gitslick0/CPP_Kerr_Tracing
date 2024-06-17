cmake_minimum_required(VERSION 3.10)

#Set the project name
project(KERR_SIM VERSION 1.0 LANGUAGES CXX Fortran)

#Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#Add the executable
add_executable(KERR_SIM ${CMAKE_SOURCE_DIR}/src/OpenGL/main.cpp ${CMAKE_SOURCE_DIR}/src/ynogk.f90)

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
target_include_directories(KERR_SIM PUBLIC ${CMAKE_SOURCE_DIR}/include)
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
message(STATUS "GLEW libraries : ${GLEW_LIBRARIES}")
message(STATUS "GLUT libraries: ${GLUT_LIBRARIES}")

set_target_properties(KERR_SIM PROPERTIES 
     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/OpenGL"
#    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build/OpenGL"
#    ARCHIVE_OUTPUT_DIRECTORY "{CMAKE_SOURCE_DIR}/build/OpenGL"
)

#-------------------------------Write Tests ------------------------------------------


# add test subdirectory
#add_subdirectory(./test)

# add testing executable
add_executable(runOpenGLTests test/main.cpp test/test_fortran_functions.cpp ${CMAKE_SOURCE_DIR}/src/ynogk.f90)
target_include_directories(runOpenGLTests PUBLIC ${CMAKE_SOURCE_DIR}/include PRIVATE ${GLFW_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS} ${OPENGL__INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} ${GLU_INCLUDE_DIRS})
target_link_libraries(runOpenGLTests gtest gtest_main glfw ${GLEW_LIBRARIES} ${OPENGL_LIBRARIES} ${GLU_LIBRARIES} ${GLUT_LIBRARIES} m)
target_compile_definitions(KERR_SIM PRIVATE
    VERTEX_SHADER_PATH="${CMAKE_SOURCE_DIR}/include/Kerr_Sim/Shaders/vShader.vs"
    FRAGMENT_SHADER_PATH="${CMAKE_SOURCE_DIR}/include/Kerr_Sim/Shaders/fShader.fs"
)

add_test(NAME OpenGLTests COMMAND runOpenGLTests)

set_target_properties(runOpenGLTests PROPERTIES 
     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/OpenGL"
#    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build/OpenGL"
#    ARCHIVE_OUTPUT_DIRECTORY "{CMAKE_SOURCE_DIR}/build/OpenGL"
)


enable_testing()