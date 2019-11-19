cmake_minimum_required (VERSION 2.6)

#>> Set up Project and Subdirectories

enable_language(Fortran)
project (ast_optim)

include_directories(${CMAKE_BINARY_DIR}/../mod/)
link_directories(${CMAKE_BINARY_DIR}/../lib/)

set (CMAKE_fortran_COMPILER gfortran)
set (CMAKE_C_COMPILER "/usr/local/gcc-8")
set (CMAKE_CXX_COMPILER "/usr/local/g++-8")
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
                #                src/jiggler.f90
                src/Fortran_Lambert_3435539_contiguous.f90
                src/state_finder_imp_jiggle.f90
                src/precision.f90
                src/constants.f90

                #src/tester.f90

)

#>> Compilation Options

target_compile_options (ast_optim PUBLIC "-Wall")
target_compile_options (ast_optim PUBLIC "-pedantic")
target_compile_options (ast_optim PUBLIC "-O3") 				                        # Optimization level -O2
target_compile_options (ast_optim PUBLIC "-ffree-line-length-512")               # Increase max line length
target_link_libraries(ast_optim PUBLIC fat spice)                               # Link the libraries
