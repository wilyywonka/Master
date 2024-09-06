cmake_minimum_required(VERSION 3.5)

project(ActiveSolidsProject VERSION 0.1
        DESCRIPTION "Active Solids Simulation Project"
        )

enable_language(Fortran)

set(CMAKE_Fortran_COMPILER gfortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/run)

add_subdirectory(source)


