# Simula
Monte Carlo Simulation for Molecular Self-Assembly (Possibly using CUDA and OpenCL for accerelation)

## How to build the program

### Linux/OSX
    cd fortran
    mkdir build
    cd build
    cmake ..
    make -j8
    ./simula

### Windows
    * build a fortran project in Visual Studio
    * add every thing inside fortran/main/src into the project
    * compile directly
