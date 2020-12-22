cmake_minimum_required (VERSION 2.6)

#>> Set up Project and Subdirectories

project (ast_optim LANGUAGES Fortran C CXX)
enable_language(Fortran C CXX)

include_directories(${CMAKE_BINARY_DIR}/../mod/) 
link_directories(${CMAKE_BINARY_DIR}/../lib/)

SET (BOOST_ROOT "~/libraries/boost")

FIND_PACKAGE(Boost ${BOOST_MIN_VERSION} REQUIRED)
find_package (Boost REQUIRED)
include_directories( ${Boost_INCLUDE_DIR} )

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

#>> Source Code

add_executable (ast_optim
                src/compute_spline.f90
                src/constants.f90
                src/integrator.cc
                src/integrator.f90
                src/main.f90
                src/midaco.o
                src/midaco_interface.f90
                src/midaco_wrap.c
                src/problem_parameters.f90
                src/sorting_routines.f90
                src/state_determination.f90
                src/variable_initialisation.f90
)

#>> Compilation Options

# Fortran
target_compile_options(ast_optim PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-Wall -O2 -cpp -DDEBUG -Wno-tabs>)
target_compile_options(ast_optim PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-O2>)
link_directories(${CMAKE_BINARY_DIR}/../lib/)
target_link_libraries(ast_optim PUBLIC fat spice bspline)  
set_property(TARGET ast_optim PROPERTY LINKER_LANGUAGE Fortran)