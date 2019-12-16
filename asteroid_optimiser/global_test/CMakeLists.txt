cmake_minimum_required (VERSION 2.6)

#>> Set up Project and Subdirectories

enable_language(Fortran)
project (ast_optim)

include_directories(${CMAKE_BINARY_DIR}/../mod/)
link_directories(${CMAKE_BINARY_DIR}/../lib/)

set (CMAKE_fortran_COMPILER "gfortran")
set (CMAKE_C_COMPILER "gcc")
set (CMAKE_CXX_COMPILER "g++")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  ENDIF()
ENDIF()

#>> Source Code

add_executable (ast_optim

                # Main project file
                src/precision.f90 
                src/state_determination.f90
                src/constants.f90
                src/main.f90
                #                src/problem.f90
                #                src/fit.f90 
                #                src/fit.f90
                #                src/testfunc.f90
                src/global.f90 
)

#>> Compilation Options

target_compile_options (ast_optim PUBLIC "-Wall")
target_compile_options (ast_optim PUBLIC "-pedantic")
target_compile_options (ast_optim PUBLIC "-O2") 				                        # Optimization level -O2
target_compile_options (ast_optim PUBLIC "-Wno-tabs")
target_link_libraries(ast_optim PUBLIC fat spice)                               # Link the libraries
