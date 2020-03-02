cmake_minimum_required (VERSION 2.6)

#>> Set up Project and Subdirectories

enable_language(Fortran)
project (ast_optim)

include_directories(${CMAKE_BINARY_DIR}/../mod/)
link_directories(${CMAKE_BINARY_DIR}/../lib/)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

#>> Source Code

add_executable (ast_optim

                # Main project file
                src/precision.f90 
                src/state_determination.f90
                src/constants.f90
                src/main.f90
                src/problem_parameters.f90
                src/compute_spline.f90
                src/variable_initialisation.f90
                src/midaco.f
                src/ancillary_data.f90
		src/mpi_variables.f90
)

#>> Compilation Options

## GNU OPTIONS
#	target_compile_options (ast_optim PUBLIC "-Wall")
#	target_compile_options (ast_optim PUBLIC "-O2") 				                        # Optimization level -O2
#	target_compile_options (ast_optim PUBLIC "-cpp")
#	target_compile_options (ast_optim PUBLIC "-DDEBUG")
#	target_compile_options (ast_optim PUBLIC "-Wno-tabs")


## IFORT OPTIONS
target_compile_options (ast_optim PUBLIC "-warn")
#target_compile_options (ast_optim PUBLIC -check all)
target_compile_options (ast_optim PUBLIC "-O2")
target_compile_options (ast_optim PUBLIC "-fpp") # Enable pre-processor
target_compile_options (ast_optim PUBLIC "-DDEBUG")
target_compile_options (ast_optim PUBLIC "-check all, noarg_temp_created")

target_link_libraries(ast_optim PUBLIC fat spice bspline-fortran)                               # Link the libraries
