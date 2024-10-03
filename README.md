# Active Solids Simulator

Active Solids Simulator is a multi-language repository consiting of high performance code used to simulate Active Solids.

### Development notes
This is still a repository that is under development, and there are multiple planned aspects and functionalities that are not yet implemented. These aspects include, but are not limited to:

- Parallel Fortran simulation code
- CUDA Fortran simulation code
- Detection of topological defects
- More analysis of the simulated systems

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

Julia is downloaded and installed by following the steps presented in this [link](https://julialang.org/downloads/platform/). The libraries used can be installed by using the Julia package manager `Pkg`. For example, if we are to install the library `Plots.jl`

```julia
import Pkg

Pkg.add("Plots")
```

After installing the libraries the Julia programs should run without issues.

### Fortran

The compillation and execution of the Fortran code is only currently tested using Linux Mint, and the following instructions may only be applicable to Debian-based OS's. The following steps should however be able to be replicated on most operating systems by modifying the steps.

To install and compile the Fortran code, there are a few prerequisites. The program uses the library [h5fortran](https://github.com/geospace-code/h5fortran) to easily read and write to HDF5-files. The installation of this is well documented on their GitHub page.

To compile and use both this code and the h5fortran library, a few HDF5 libraries have to be installed and linked to as well, and the paths to these libraries needs to be set up. In the file located at `FortranXXXX/source/CMakeListsTemplate.txt` there are multiple places where the placeholder text `!FILL!` is written. Here the path to the relevant library must be written in place of the placeholder `!FILL!`, before renaming the file to `FortranXXXX/source/CMakeLists.txt`, where `FortranXXXX` indicates the three folders in the repository, Serial, Parallel and CUDA.
(The template was made such that the machine specific paths were not set, as to avoid confusion.)

After correcting these paths the program can be compiled. CMake is chosen as the compilation tool and the program can be compiled by executing

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

The parameters read by the program will then be written to the screen, and the result will be written to the HDF5-file with the specified savefile name.

### A note on HDF5-files
*Tl;dr: If the same file name is multiple times to save a HDF5 file, like `SaveFiles.h5`, delete the existing file before a new simulation.*

When choosing a filename for the file containing the snapshots, the program does not delete existing files with the same name. It will instead partially write over this file, due to the nature of HDF5-files. 
This can lead to two main problems, with the first being that the program can crash if e.g. the number of paricles change between runs, and the file is still present. The reason is that the program will attempt to write over an array that already exists, but as the dimensions now have changed, the operation will fail.
The second can occur when only the number of snapshots is reduced between two runs. The program will run without errors, and the snapchots will be saved to the savefile. The problem is that as it is now overwriting existing arrays, the previous simulation snapshots will still exist. This will lead to incorrect results if the analysis attempts to analyze both the current snapshots, and some left over snapshots from a previous run.


## Creator and acknowledgements

This software was created by Wilhelm Sunde Lie as a part of the Master's Thesis in Physics at NTNU

My thanks go to my supervisor Prof. Paul Gunnar Dommersnes for a lot of the physics related questions, and to my co-supervisor Prof. Ingve Simonsen for a lot of the programming/Fortran/compilation questions.

## License
This project is licensed under the [MIT License](https://choosealicense.com/licenses/mit/).