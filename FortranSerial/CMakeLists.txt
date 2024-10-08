cmake_minimum_required(VERSION 3.5)

project(ActiveSolidsProject VERSION 0.1
        DESCRIPTION "Active Solids Simulation Project"
        )

enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -xHOST")# -Ofast -xHOST
    set(CMAKE_Fortran_FLAGS_PROFILE "-Ofast -xHOST -g")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)    
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall") # -Wall
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native")# -Ofast -march=native
    set(CMAKE_Fortran_FLAGS_PROFILE "-pg")
endif()


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/run)

add_subdirectory(source)

 
