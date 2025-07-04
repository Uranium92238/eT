# Prerequisites
1. [CMake](https://cmake.org/) (3.7 or newer)
2. Python 3 (3.6 or newer)
3. Recent GNU (gfortran, gcc, g++ 8 or newer) or Intel compilers (ifort*, icc, icpc 2020 or newer)
4. BLAS and LAPACK libraries
5. [Libint 2 library](https://github.com/evaleev/libint)
with integrals for one-body operators and electron repulsion enabled.
Libint has the dependencies [Eigen 3](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [Boost](https://www.boost.org).
To compile Libint, follow the instructions below or consult the [Libint Wiki](https://github.com/evaleev/libint/wiki).


## Installing Libint
Download the [Libint library for eT](https://www.etprogram.org/libint/libint-2.7.0-beta.6.tgz).
Unpack the tar file:
```shell
tar -xvzf libint-2.7.0-beta.6.tgz
```
Enter the generated folder:
```shell
cd libint-2.7.0-beta.6
```
Compile (with 4 threads (-j flag), increase if more are available):
```shell
cmake . -DCMAKE_INSTALL_PREFIX=/where/you/want/to/install/libint/libint-2.7.0-beta.6 -DCMAKE_C_COMPILER=[C compiler] -DCMAKE_CXX_COMPILER=[C++ compiler] -DCMAKE_CXX_FLAGS=[C++ compiler flags]
cmake --build . -j 4
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
Also note that `--prefix=/where/you/want/to/install/pcmsolver/` should be used
as command line flag when running `setup.py` for PCMSolver.
CMake should then be able to locate it,
if the location of PCMSolver is set as the following environment variable:
```shell
export PCMSolver_ROOT=/where/you/want/to/install/pcmsolver/
```
On the [wiki](https://gitlab.com/eT-program/eT/-/wikis/home),
we have a
[bash script](https://gitlab.com/eT-program/eT/-/wikis/Various%20guides/Example%20script%20for%20installing%20PCMSolver)
that can be used to setup PCMSolver.



# Installing eT
## Getting the eT source
The source code can be downloaded as a tarball from the
[releases](https://gitlab.com/eT-program/eT/-/releases) page.
Click the  `source code` dropdown menu and choose your favourite archive format.
After downloading,
your archive manager will probably offer to handle it.
You can also use the command line.
```shell
unzip eT-vx.y.z.zip
tar -xvzf eT-vx.y.z.tar.gz
tar -xvjf eT-vx.y.z.tar.bz2
tar -xvf eT-vx.y.z.tar
```
Note that the archives do not include the
[runtest](https://runtest.readthedocs.io/en/latest/)
submodule required for testing
or the optional
[pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit)
submodule to run unit tests.

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
It will automatically download the submodules
[runtest](https://runtest.readthedocs.io/en/latest/),
a python program used to run the test suite,
and
[pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit),
to run unit tests.
Note that pFUnit required CMake version 3.12 or newer.
If you don't use `--recursive`,
but later change you mind,
you can get the submodule with the following commands.
```shell
git submodule init
git submodule update
```

## Configuration of eT
After downloading,
you will have a directory called `eT` or something like `eT-v1.0.1`,
depending on how you did it.
We will refer to it as `eT` from now on.
Go to this directory:
```shell
cd eT
```
Run the [setup.py](https://etprogram.org/setup.html) script to configure CMake:
```shell
./setup.py
```
Compilers and libraries identified by CMake will be printed to screen.
Take a look and see if it looks reasonable and error free.
Note that `setup.py` may take a long time on some clusters.

CMake will try to identify the compilers to use (Fortran, C, C++) and the location of libraries (BLAS, LAPACK, Libint, Eigen, Boost).
If CMake does not correctly locate and identify these,
try setting the associated environment variables.
In order to have these variables automatically exported when you open a new shell,
you can place the export commands in your `.bashrc` file.
Remember to update your current shell with `source .bashrc`.
Below are some examples:
```shell
export LIBINT2_ROOT=/home/username/prog/libint
export Eigen3_ROOT=/home/username/prog/eigen3/include/eigen3
export BOOST_INCLUDEDIR=/usr/include
export MATH_ROOT=/opt/intel/mkl
```
If CMake still does not find Eigen3,
it can help to add `EIGEN3_INCLUDE_DIR` to the `.bashrc`
and also add it to the cmake prefix path:
```shell
export EIGEN3_INCLUDE_DIR=/usr/local/include/
export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$EIGEN3_INCLUDE_DIR"
```
Similarly, Boost can be added to the prefix path as well, if it cannot be found
```shell
export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$BOOST_INCLUDEDIR"
```

For help with `setup.py`,
run the script with the `--help` option.
```shell
./setup.py --help
```
This will list the various options,
a short description,
and the default values.
For more detailed description,
see the [website](https://etprogram.org/setup.html).

**Optional:**
To enable PCMSolver,
you must run `setup.py` with the `--pcm` option.

## Compilation of eT
If everything went well,
you should now have a directory in eT called build.
Go to this directory and use the `make` command to compile.
```shell
cd build
make
```
Note that either the `setup.py` or `autogenerate_files.py` script have to be run before
compiling, as the script generates necessary interfaces and files handling complex variables in eT.
Use `make -j n`, where `n` is the number of processors, to make the compilation run in parallel.

If successful, the directory now contains the executable, `eT`, as well as a Python launch script
[eT_launch.py](https://etprogram.org/eT_launch.html).
If you downloaded runtest,
you can test your installation using ctest:
```shell
ctest
```

## Running eT
The launch script is the recommended way of running eT.
Similarly to `setup.py`,
it has a `--help` option
and a more detailed description is available
[here](https://etprogram.org/eT_launch.html).
To have `eT_launch.py` easily available,
you can add `build` to your `PATH` variable in `.bashrc`;
```shell
export PATH=$PATH:/path/to/eT/build
```
or copy the script to somewhere more convenient.

See the [website](https://etprogram.org)
for more help or look in `eT/tests` for inspiration when making your own input file.
For example,
`eT/tests/hf_energy/hf_energy.inp` is an input file for a Hartree-Fock calculation.
