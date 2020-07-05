cmake_minimum_required (VERSION 2.6)

#>> Set up Project and Subdirectories

enable_language(Fortran)
project (pre_filter)

link_directories(${CMAKE_BINARY_DIR}/../lib/)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

#>> Source Code

add_executable (pre_filter

                # Main project file
                src/precision.f90 
                src/constants.f90
                src/main.f90
                src/problem_parameters.f90
                src/variable_initialisation.f90
                src/utility_functions.f90
)

#>> Compilation Options

target_compile_options (pre_filter PUBLIC "-Wall")
target_compile_options (pre_filter PUBLIC "-O2")
target_compile_options (pre_filter PUBLIC "-g")
target_compile_options (pre_filter PUBLIC "-cpp")
target_compile_options (pre_filter PUBLIC "-Wno-tabs")
target_compile_options (pre_filter PUBLIC "-Wno-array-temporaries")
target_link_libraries  (pre_filter PUBLIC spice)
