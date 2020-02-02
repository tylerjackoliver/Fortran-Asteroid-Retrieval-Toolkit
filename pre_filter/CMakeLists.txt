cmake_minimum_required (VERSION 2.6)

#>> Set up Project and Subdirectories

enable_language(Fortran)
project (pre_filter)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

#>> Source Code

add_executable (pre_filter

                # Main project file
                src/precision.f90 
                src/constants.f90
                src/main.f90
                src/problem_parameters.f90
                src/variable_init.f90
)

#>> Compilation Options

target_compile_options (pre_filter PUBLIC "-Wall")
target_compile_options (pre_filter PUBLIC "-O2") 				                        # Optimization level -O2
target_compile_options (pre_filter PUBLIC "-Wno-tabs")
