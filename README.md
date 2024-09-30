# Active Solids Simulator

Active solids simulator is a multi-language repository consiting of high performance software used to simulate active solids.

### Development notes
This is still a repository that is under development, and there are multiple planned aspects and functionalities that are not yet implemented. These aspects include, but are not limited to:

- Parallel Fortran simulation code
- CUDA Fortran simulation code
- Detection of topological defects

The author achknowledges that there may very well be mistakes in the code during development and testing.

## Description

The repository consists of two main languages, Julia and Fortran.

The Julia program `InitSetup.jl` is used for creating the initial systems by simulating the growth of differently sized discs, and writing this initial system to a HDF5-file, as well as creating a parameter NameList to be read by the Fortran program.

The Fortran program `ActiveSolidsSimulation` simulates the active solid using parameters set in the NameList `Parameters.nml` and the initial system saved in `InitState.h5`. It will then save snapshots of the system to the HDF5 file `SaveFiles.h5`.

After Fortran program is run, the Julia program `Analysis.jl` is run to analyze and visualize the simulated system.

## Installation

The project is divided into four separate folders: Julia, FortranSerial, FortranParallel and FortranCUDA.

Start by cloning the repository to your local machine.

### Julia 


### Fortran

Before 


## Creator and acknowledgements

This software was created by Wilhelm Sunde Lie as a part of the Master's Thesis in Physics at NTNU

My thanks goes to to my supervisor Prof. Paul Gunnar Dommersnes for a lot of the physics related questions, and to my co-supervisor Prof. Ingve Simonsen for a lot of the programming/Fortran/compilation questions.

## 