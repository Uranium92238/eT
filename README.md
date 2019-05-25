# Prerequisites
1. [Libint 2.4.2 library](https://github.com/evaleev/libint/releases/download/v2.4.2/libint-2.4.2.tgz)
2. [Eigen 3](http://eigen.tuxfamily.org/index.php?title=Main_Page) (3.3.5 or newer)
3. [CMake](https://cmake.org/) (3.10 or newer)
4. GNU compilers (gfortran, gcc, g++) or Intel compilers (ifort, icc, icpc). These need to be recent versions (2016-) since some features of the Fortran 2008 standard were included only recently in Gfortran and Ifort. 
5. BLAS and LAPACK libraries.

You also need to set two environment variables, e.g. in .bashrc:
* export LIBINT\_DATA\_PATH=/path-to-libint/libint-2.4.2/lib/basis
* export SAD\_ET\_DIR=/path-to-eT/src/molecular\_system/sad

# Installation of eT
Open the terminal and clone the repository:
```
git clone git@gitlab.com:eT-program/eT.git
```
This will create a folder called "eT". Change to this directory:
```
cd eT
```
Initialize and update the submodules:
```
git submodule init 
git submodule update
```
Now make a directory for the compiled program, and enter the directory:
```
mkdir build
cd build
```
Prepare the files with CMake:
```
cmake ..
```
Compile eT:
```
make
```
Run tests:
```
ctest
```
