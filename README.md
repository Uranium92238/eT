# Prerequisites
1. [CMake](https://cmake.org/) (3.7 or newer)
2. Python 3 (3.5 or newer)
3. Recent (2016-) GNU (gfortran, gcc, g++) or Intel compilers (ifort, icc, icpc) 
4. BLAS and LAPACK libraries
5. [Libint 2 library](https://github.com/evaleev/libint) 
with integrals for one-body operators and electron repulsion enabled.
Libint has the dependencies [Eigen 3](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [Boost](https://www.boost.org). 
To compile Libint, follow the instructions below or consult the [Libint Wiki](https://github.com/evaleev/libint/wiki),
or see below.

## Installing Libint
Download the [Libint library for eT](https://www.etprogram.org/libint/libint-2.7.0-beta.1.tgz). 
Unpack the tar file:
```shell
tar -xvzf libint-2.7.0-beta.1.tgz
```
Enter the generated folder:
```shell
cd libint-2.7.0-beta.1
```
Compile:
```shell
cmake . -DCMAKE_INSTALL_PREFIX=/where/you/want/to/install/libint/libint-2.7.0-beta.1 -DCMAKE_C_COMPILER=[C compiler] -DCMAKE_CXX_COMPILER=[C++ compiler] -DCMAKE_CXX_FLAGS=[C++ compiler flags]
cmake --build .
```
CMake will attempt to install Libint in the directory specified by `-DCMAKE_INSTALL_PREFIX`. 
Especially on clusters this is important, 
since you normally won't have access to the default location, `/usr/local`. 
Provide a path to a directory you want to install in and make sure you have write access.
The install directory cannot be in the source directory.

Note that when compiling with Intel, one should include the following flags:
```shell
-DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_CXX_FLAGS=-std=c++11
```
Make sure the Libint configure step finds the correct compilers (as compiling Libint can take several hours). Compile and install:
```shell
cmake --build . --target install
```

On a cluster, compiling Libint will often take so long that you must (and want to) submit it as a job.
On the [wiki](https://gitlab.com/eT-program/eT/-/wikis/home), 
we have an [example script](https://gitlab.com/eT-program/eT/-/wikis/Various-guides/Example-script-for-installing-Libint-on-cluster) to submit the compilation. 
You will probably have to modify it for your system.
 
If you wish to download and run the Libint compiler yourself, 
the standard configuration flags for eT are available [here](https://gitlab.com/eT-program/eT/-/wikis/Various-guides/Standard-configuration-for-Libint).


# Optional: Installing PCMSolver
If you want to use polarizable continuum models in eT, 
you can link eT to the [PCMSolver](https://github.com/PCMSolver/pcmsolver) library.
PCMSolver can be cloned from Github and must be compiled separately. 
Instructions for downloading and compiling can be found [here](https://pcmsolver.readthedocs.io/en/stable/). 
A version that compiles on Mac and passes all eT tests with array bound checks is available [here](https://github.com/eirik-kjonstad/pcmsolver).
Note that PCMSolver requires [Zlib](https://www.zlib.net/) and [Boost](https://www.boost.org).



# Installing eT
## Getting the eT source
The source code can be downloaded as a tarball from the 
[releases](https://gitlab.com/eT-program/eT/-/releases) page. 
Click the  `source code` dropdown menu and choose your favourite archive format.
After downloading,
your archive manager will probably offer to handle it. 
You can also use the command line.
```shell
unzip eT-v1.0.0.zip
tar -xvzf eT-v1.0.0.tar.gz
tar -xvjf eT-v1.0.0.tar.bz2
tar -xvf eT-v1.0.0.tar
```
Note that the archives do not include the runtest submodule required for testing.

You can also download eT by opening a terminal and cloning the repository:
```shell
git clone --recursive https://gitlab.com/eT-program/eT.git 
```
or
```shell
git clone --recursive git@gitlab.com:eT-program/eT.git 
```
The second option requires a user account on [Gitlab](https://gitlab.com/),
and that you have set up a public/private [keypair](https://docs.gitlab.com/ee/ssh/),
but is the most convenient if you intend to contribute.

**Note:**
`--recursive` is optional, but recommended. 
It will automatically download [runtest](https://runtest.readthedocs.io/en/latest/), 
a python program used to run the test suite.
If you don't use `--recursive`, 
but later change you mind, 
you can get the submodule with the following commands.
```shell
git submodule init
git submodule update
```

## Configuration of eT
After downloading, 
you will have a directory called `eT` or `eT-v1.0.0`,
depending on how you did it. 
We will refer to it as `eT` from now on.
Go to this directory:
```shell
cd eT
```
Run the [setup](https://gitlab.com/eT-program/eT/-/wikis/Using-eT/How-to-use-the-setup-script) script to configure CMake:
```shell
./setup 
```
Compilers and libraries identified by CMake will be printed to screen.
Take a look and see if it looks reasonable and error free.
Note that `setup` may take a long time on some clusters.

CMake will try to identify the compilers to use (Fortran, C, C++) and the location of libraries (BLAS, LAPACK, Libint, Eigen, Boost). 
If CMake does not correctly locate and identify these, 
try setting the associated environment variables.
In order to have these variables automatically exported when you open a new shell, 
you can place the export commands in your `home/.bashrc` file. 
Remember to update your current shell with `source .bashrc` Below are some examples:
```shell
export LIBINT2_ROOT=/home/username/prog/libint
export Eigen3_ROOT=/home/username/prog/eigen3/include/eigen3
export BOOST_INCLUDEDIR=/usr/include
export MATH_ROOT=/opt/intel/mkl
```

For help with `setup`, 
run the script with the `--help` option.
```shell
./setup --help
```
This will list the various options, 
a short description,
and the default values.
For more detailed description, 
see the [wiki page](https://gitlab.com/eT-program/eT/-/wikis/Using-eT/How-to-use-the-setup-script).

**Optional:**
To enable PCMSolver, 
you must run `setup` with the `--pcm` option.
It might be a good idea to set the location of PCMSolver as an environment variable to help CMake locate it.
```shell
export PCMSolver_ROOT=/path/to/pcmsolver
```

## Compilation of eT
If everything went well,
you should now have a directory in eT called build.
Go to this directory and use the `make` command to compile.
```shell
cd build
make
```
Use `make -j n`, where `n` is the number of processors, to make the compilation run in parallel.

If successful, the directory now contains the executable, `eT`, as well as a Python launch script 
[eT_launch](https://gitlab.com/eT-program/eT/-/wikis/Using-eT/How-to-use-the-launch-script). 
If you downloaded runtest,
you can test your installation using ctest:
```shell
ctest
``` 

## Running eT
The launch script is the recommended way of running eT. 
Similarly to `setup`, 
it has a `--help` option 
and a more detailed description is available 
[here](https://gitlab.com/eT-program/eT/-/wikis/Using-eT/How-to-use-the-launch-script).
To have `eT_launch` easily available,
you can add `build` to your `PATH` variable in `.bashrc`;
```shell
export PATH=$PATH:/path/to/eT/build
```
or copy the script to somewhere more convenient.

See the [wiki](https://gitlab.com/eT-program/eT/-/wikis/home) 
for more help or look in `eT/tests` for inspiration when making your own input file.
For example, 
`eT/tests/hf_energy/hf_energy.inp` is an input file for a Hartree-Fock calculation.

