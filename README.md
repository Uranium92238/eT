# Prerequisites
1. [CMake](https://cmake.org/) (3.7 or newer)
2. Python 3
3. Recent (2016-) GNU (gfortran, gcc, g++) or Intel compilers (ifort, icc, icpc) 
4. BLAS and LAPACK libraries
5. [Libint 2 library](https://www.etprogram.org/libint/libint-2.7.0-beta.1.tgz) (link: 2.7.0 library configured for eT)
and the dependencies [Eigen 3](http://eigen.tuxfamily.org/index.php?title=Main_Page) (3.3.5 or newer) and [Boost](https://www.boost.org). 
6. To compile the library use the following commands, or see the instructions on the [Libint Wiki](https://github.com/evaleev/libint/wiki).

## Install Libint
Download the Libint library for eT. Unpack the library:
```
tar -xvzf libint-2.7.0-beta.1.tgz
```
Enter the generated folder:
```
cd libint-2.7.0-beta.1
```
Compile:
```
cmake . -DCMAKE_INSTALL_PREFIX=/usr/local/libint/libint-2.7.0-beta.1 -DCMAKE_CXX_COMPILER=[C++ compiler] CXXFLAGS=[C++ compiler flags]
cmake --build .
```
In the compilation step, the installation prefix should be provided.


Install:
```
cmake --build . --target install
```
 
# Quick install of eT
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
Run the setup script to configure using CMake:
```
./setup [--help]
```
CMake will identify the compilers to use (Fortran, C++) and the location of libraries (BLAS, LAPACK, Libint, Eigen, Boost). If the setup script does not correctly locate and identify libraries (BLAS, LAPACK, Libint, Eigen, Boost), try setting the associated environment variables in the .bashrc (example paths given):
```
export LIBINT2_ROOT=/home/eirikfad/prog/libint
export EIGEN3_ROOT=/home/eirikfad/prog/eigen3/include/eigen3
export BOOST_INCLUDEDIR=/usr/include
export MATH_ROOT=/opt/intel/mkl/lib/intel64_lin 
```
Change to the "build" directory and compile the executable:
```
cd build
make [-j4]
```
If successful, the directory now contain the executable (eT) as well as a Python launch script (eT_launch). 
To test that the program performs as expected, make sure all the tests pass using the compiled executable:
```
ctest
``` 
To use the launch script, it is useful to define an alias in .bashrc:
```
alias eT=pathtoeT/build/eT_launch
```
The program can then be conveniently run by using the command eT:
```
eT [--help]
```
