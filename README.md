# Quick install of eT

Open the terminal and clone the repository:
```
git clone git@gitlab.com:eT-program/eT.git
```
This will create a folder called "eT". Change to this directory:
```
cd eT
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
Note: eT is written in Fortran 2008 standard and requires very recent versions of CMake and Gfortran (or Ifort) compilers that supports the use of submodules.

<< how to install libint here >>
