cmake_minimum_required(VERSION 3.10)

#Set the project name
project(FITS_WRITE VERSION 1.0 LANGUAGES CXX Fortran)

#Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#Add the executable
add_executable(FITS_WRITE ${CMAKE_SOURCE_DIR}/src/FITS_WRITE/main.cpp)

find_package(CFITSIO REQUIRED)

include_directories(${CFITSIO_INCLUDE_DIR})

target_link_libraries(FITS_WRITE ${CFITSIO_LIBRARIES})

set_target_properties(FITS_WRITE PROPERTIES 
     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/FITS_WRITE"
#    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build/OpenGL"
#    ARCHIVE_OUTPUT_DIRECTORY "{CMAKE_SOURCE_DIR}/build/OpenGL"
)

add_executable(runFITSTests test/main.cpp test/Fitstest.cpp)
target_include_directories(runFITSTests PUBLIC ${CMAKE_SOURCE_DIR}/include PRIVATE ${CFITSIO_INCLUDE_DIR})
target_link_libraries(runFITSTests gtest gtest_main ${CFITSIO_LIBRARIES})

add_test(NAME FITSTests COMMAND runFITSTests)

set_target_properties(runFITSTests PROPERTIES 
     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/FITS_WRITE"
#    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build/OpenGL"
#    ARCHIVE_OUTPUT_DIRECTORY "{CMAKE_SOURCE_DIR}/build/OpenGL"
)


enable_testing()