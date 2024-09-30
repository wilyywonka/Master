# Active Solids Simulator

Active solids simulator is a multi-language repository consiting of high performance code used to simulate Active Solids.

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

## Installation and Use

The project is divided into four separate folders: Julia, FortranSerial, FortranParallel and FortranCUDA.

Start by cloning the repository to your local machine.

### Julia 

Julia is downloaded and installed by following the steps presented in this \[link]{https://julialang.org/downloads/platform/}. The libraries used can be installed by using the julia package manager `Pkg`. For example, if we are to install the library `Plots.jl`

```julia
import Pkg

Pkg.add("Plots")
```


### Fortran

The program is only currently tested using Linux Mint.

To install and compile the Fortran code, there are a few prerequisites. The program uses the package \[h5fortran]{https://github.com/geospace-code/h5fortran} to easily read and write to HDF5-files. The installation of this is well documented on their GitHub page.

To compile and use both this code and the h5fortran library, a few HDF5 libraries have to be installed and linked to as well, and the paths to these libraries are set in the `FortranXXXX/source/CMakeLists.txt`, where `FortranXXXX` indicates the three folders in the repository, Serial, Parallel and CUDA.

After setting these paths the program can be compiled. CMake is chosen as the compilation tool and the program can be compiled by executing

```bash
cd FortranXXXX/build

cmake ..

make
```

where again `FortranXXXX` is replaced by either of the three folder names. The program will be located in the `FortranXXXX/build/run` and can be run by executing

```bash
cd run

./ActiveSolidsSimulation
```

The parameters read by the program will then be written to the screen.

*** NOTE: *** Tl;dr: Delete the `SaveFiles.h5` before a new simulation. The Program does not delete the previous `SaveFiles.h5` file automatically, and can fail if e.g. the number of paricles change between runs, and the file is still present. Due to the nature of HDF5-files, if only the number of snapshots is reduced after a run, the program will run without errors, and the snapchots will be saved to the `SaveFiles.h5` file, but previous simulation snapshots will still exist in the file, and will not be written over. This can lead to weird results if more than one simulation's snapshots are analyzed.


## Creator and acknowledgements

This software was created by Wilhelm Sunde Lie as a part of the Master's Thesis in Physics at NTNU

My thanks goes to to my supervisor Prof. Paul Gunnar Dommersnes for a lot of the physics related questions, and to my co-supervisor Prof. Ingve Simonsen for a lot of the programming/Fortran/compilation questions.

## 