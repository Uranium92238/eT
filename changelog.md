# eT v1.4.0
### Bugfixes
- Fixed some errors in memory handling that can occur for large systems and/or on machines with limited memory resources. eT-program/eT!852

### Features
- Timers for Fock matrix construction for coupled cluster. eT-program/eT!849
- `cc_propagation` and `complex_fft` file names are now starting with `output_name.*` and automatically copied to the output directory. eT-program/eT!854
- Added option to plot CNTOs and NTOs. eT-program/eT!816
- Runtime check for reading of non-existing keywords/sections. eT-program/eT!872
- Runtime check that keywords have been defined for sections. eT-program/eT!878
- TDHF (RPA/Tamm-Dancoff) excitation energies implemented for RHF. eT-program/eT!809
- Restart for TDHF (RPA/Tamm-Dancoff) excitation energies. eT-program/eT!883
- Added support for unit testing with pFUnit package. eT-program/eT!787
- Minor changes to memory-handling to improve testability of memory-batching. eT-program/eT!852 
- New layout for timings.out file eT-program/eT!891

### Optimization
- Adding option for overlap screening to CC calculation of one electron integral. eT-program/eT!851
- Record storers now delete their file when finalized is called. eT-program/eT!893

### Structure
- Added block class using the range class. Facilitates the handling of multiple index ranges. eT-program/eT!823
- Angular momentum tools handle components of the angular momenta. eT-program/eT!864

# eT v1.3.9
### Bugfixes
- Updated Readme with `setup.py` and `eT_launch.py`. eT-program/eT!887

# eT v1.3.8
### Bugfixes
- Fix bug with pure gaussians. eT-program/eT!884

# eT v1.3.7
### Bugfixes
- Remove convergence check for cc2 test with davidson preconvergence. eT-program/eT!881

# eT v1.3.6
### Bugfixes
- Basis specification for atom subsets is now possible again. eT-program/eT!880

# eT v1.3.5
### Bugfixes
- Missing keyword (one-electron integral cutoff) added to the input file class. eT-program/eT!874

# eT v1.3.4
### Bugfixes
- Bug with reordering of atoms with both active spaces and several different basis sets was fixed. eT-program/eT!868

# eT v1.3.3
### Bugfixes
- Ordering of AOs for spherical gaussians fixed in the orbital files. eT-program/eT!857
- Molden file now contains correctly normalized orbital coefficients. eT-program/eT!857

# eT v1.3.2
### Bugfixes
- Read response vector tk for polarizabilities when needed for either diagonal or cross terms, to avoid cross terms from uninitialized array. eT-program/eT!856

# eT v1.3.1
### Bugfixes
- Newlines added in release script for correct formatting of release notes. eT-program/eT!846

# eT v1.3.0
### Bugfixes
- Removed option to accumulate and erase history in DIIS. These features were not working properly. eT-program/eT!800
- Keyword `pure gaussians` is now correctly read and enabled. eT-program/eT!810
- Added check for 'end geometry' in input file. eT-program/eT!825
- Added check for no calculation specified on input. eT-program/eT!831

### Features
- Added multimodel Newton-Raphson capability to ground state Newton-Raphson solver. Can be used to converge amplitudes in CC2 and CC3. eT-program/eT!759
- EOM polarizabilities are now available for CC2. eT-program/eT!761
- Consistent copying of all output files. They will now be copied to the same directory as output file. eT-program/eT!763
- Possible to specify output file, executable, savedir, loaddir, loadv1.0dir and pcm input file for each input file. eT-program/eT!763
- `eT_launch.py` now more tolerant of number of inputs and will use defaults or give a warning if they don't match the number of input files. eT-program/eT!763
- Changed `convert_v1_0_restart_files.py` script to not delete the old files and it can now be called from `eT_launch.py` using `--load-v1-0-restart-dir`. eT-program/eT!763
- `-i-err` changed to  `-i`. (`-i-err` retained for backwardness). eT-program/eT!763
- Possible to tell `eT_launch.py` to gracefully terminate eT using signals (`--signals`) or timeout (`--timeout`). eT-program/eT!763, eT-program/eT!820
- Z-matrix tool can convert between xyz coordinates and Z-matrix. eT-program/eT!619
- Davidson excited state solver can now handle cases where a requested root forms a complex pair with a higher root before convergence is reached (i.e., a false-positive complex pair); we thus avoid the commonly encountered "add one more root" error. eT-program/eT!772
- Biorthonormalization works for a set of close states as well. eT-program/eT!775
- Option to write a molden file containing HF data. eT-program/eT!788
- Added printing of implementation references, and the associated keyword `full references` in the `print` section. eT-program/eT!773
- Made autogenerate script compatible with Python 3.6 (UTF-8 encoding error). eT-program/eT!794
- Added `ghost` atoms. eT-program/eT!795
- Added `remove core` projection to CC3. eT-program/eT!804
- Excited state properties and oscillator strengths between excited states can be calculated. eT-program/eT!785
- Added timings for each iteration in the excited state solvers. eT-program/eT!834
- Renamed input files in tests the directories to agree with the name of the directory. eT-program/eT!833

### Optimization
- One-electron integrals with screening based on the overlap matrix. eT-program/eT!770

### Structure
- The test scripts now uses pathlib instead of os. eT-program/eT!756
- Changed names of `setup` and `eT_launch` to `setup.py` and `eT_launch.py` so they can be imported into python. eT-program/eT!763
- Changed output file names to all start with `eT.*` and `eT_launch.py` renames them to `output_name.*`. eT-program/eT!763
- Removed `check_alloc.py`. eT-program/eT!763
- Cholesky orbital tool now handles decomposition of densities to get orbitals. eT-program/eT!762
- Fixed warnings from intel compilers, renaming `autogenerated_complex_files` to `generated_complex_files`. eT-program/eT!786
- Replaced interval class with more general range class. eT-program/eT!817
- Updated Fortran version to 2018 standard. eT-program/eT!836

# eT v1.2.5
### Bugfixes
- The `req_single_batch` in `calculate_energy_ccs` is now passed as `req_single_batch` instead of `element_size`. eT-program/eT!835

# eT v1.2.4
### Bugfixes
- Engines are passed by reference in oei.cpp to avoid copy of engine for each call to the `construct_oei` routine (done for each shell pair).eT-program/eT!821
- Fixed memory required for a single batch in `omega_ccsd_a2`. eT-program/eT!818

# eT v1.2.3
### Bugfixes
- Reduced threshold in `xenon` test to 1.0d-9 due to instabilities. See merge request eT-program/eT!782

# eT v1.2.2
### Bugfixes
- Testing that the number of electrons is even for restricted closed calculations. See merge request eT-program/eT!777

# eT v1.2.1
### Bugfixes
- Reduced threshold in `uhf_energy_cs` to 1.0d-10 due to instabilities. See merge request eT-program/eT!766

# eT v1.2.0
### Bugfixes
- Added handling of zero length character prints. See merge request eT-program/eT!716

### Features
- printf accepts complex numbers and possible to repeat formats. See merge request eT-program/eT!674
- UHF spin contamination is computed and printed in all UHF calculations. See merge request eT-program/eT!675
- Tidied up Omega and Jacobian timings to make them more consistent. See merge request eT-program/eT!684
- SAD guess can be used for elements up to Dubnium. See merge request eT-program/eT!707
- Frozen core can be used for elements up to Dubnium. See merge request eT-program/eT!709
- `eT_launch` now returns the lanczos files. See merge request eT-program/eT!710
- Gaussian .cube format for density and orbital visualization. See merge request eT-program/eT!603
- Specific cores can now be removed (by projection) in CC valence excited states. See merge request eT-program/eT!713
- Default threshold to check for degenerate states reduced to 1.0d-3. See merge request eT-program/eT!724
- Parallel excited states are now removed from the wavefunction if found eT-program/eT!717
- Added keyword to control storage in micro iterations of Newton-Raphson ground state solver (`micro iterations storage`). See merge request eT-program/eT!651
- Non-linear davidson: minimal threshold in the microiterations set to half the residual threshold. See merge request eT-program/eT!738
- The conversion script removes `scf_restart_file` and `cc_restart_file` after it has been used. See merge request eT-program/eT!745

### Optimization
- Fock matrix construction at CC level, in terms of T1-transformed integrals, is now N^4 scaling. Only the necessary blocks for a given task are constructed. See merge request eT-program/eT!678
- Removed one construction of the Coulomb fock matrix per iteration. See merge request eT-program/eT!677
- CCS Jacobian is now N^4 scaling (Cholesky vector algorithm) See merge request eT-program/eT!592
- CCS energy is now N^4 scaling (Cholesky vector algorithm) See merge request eT-program/eT!715
- CC3 optimization using covariant/contravariant amplitudes and integrals are resorted on the fly. See merge request eT-program/eT!685
- Non-linear davidson now default solver for CC3 and low-mem CC2. See merge request eT-program/eT!724
- Only necessary (vv|vv) integrals computed in CCSD and transpose Jacobian transformation now uses the optimized v^4o^2 routine. See merge request eT-program/eT!711
- Some cleanup in C2 and D2 terms in CCSD. See merge request eT-program/eT!727
- ERI tools use xsyrk to compute symmetric integrals. See merge request eT-program/eT!730

### Structure
- Omega and Jacobian routines have separate in and out vectors and call the corresponding parent routines. See merge request eT-program/eT!684
- An AO tool is added that handles the basis set, AO integral construction, and  all communication with Libint; replaces the molecular system class
- The MM and PCM code is restructured in a hierarchy of embedding/environment classes. Wave function classes and embedding classes are more decoupled.
- Input file routines for keywords and sections have been shortened, with some additional cleanup of input variables/routines. See merge request eT-program/eT!705.
- Autogenerated directories and interfaces are now ignored by git. They are generated by the `setup` or the `autogenerate_files` script. See merge request eT-program/eT!706
- A single SCF solver now handles all SCF calculations, it can be used with different types of convergence acceleration (none, DIIS, CROP). See merge request eT-program/eT!654
- Removed some unused routines from reordering and array utilities. See merge request eT-program/eT!725
- Newton-Raphson solver for coupled cluster ground state is now handled by the same solver as the default DIIS solver. See merge request eT-program/eT!651
- Optimized the autogenerate script and added functionality to generate complex modules from real modules. See merge request eT-program/eT!735


# eT v1.1.3
### Structure
- Made Python script to update all test refences and updated the files. See merge request eT-program/eT!746


# eT v1.1.2
### Bugfixes
- Fixed incorrect memory estimate in Jacobian transpose ccsd. See merge request eT-program/eT!718


# eT v1.1.1
### Features
- Functionality added for skipping SCF if equations are solved.
  Restart between 1.0.x and 1.1 now possible with the conversion script. See merge request eT-program/eT!683
- Global restart keyword for `scf geoopt`. See merge request eT-program/eT!693

### Structure
- `mo_information_file` now part of `hf_class` as workaround for ifort segfault. See merge request eT-program/eT!681


# eT v1.1.0
### Bugfixes
- Fixed bug in EOM when $`n_o > n_v`$
- Fixed wrong allocation in expectation value engine
- Fixed wrong memory estimate for doubles A1 term, CCSD A2 term, and CCSD E1 (Jacobian transpose)
- File needed for restart in MLHF no longer deleted
- Threshold in biorthonormalization no longer squared
- Fixed the linear dependence threshold for all davidson solvers, the threshold is now 1.0d-11 or lower.
- Fixed wrong dimension of optwork in eigen\_davidson\_tool to comply with gfortran10
- If $`e^T`$ is compiled with 32bit integers, it will check that the number of MOs does not exceed 215
  so that $`n_{MO}^4`$ can be stored in an integer.
- Fixed the useage of 32bit integers in the batching setup.
- Fixed a bug with intel compilers where we could not use printf to print two chars returned by two different functions.
- Fixed bug where the linear dependence threshold could be larger than the residual threshold.
- The label in the print of the largest amplitudes was changed from 'r' to 'tbar' for the multipliers.
- Fixed basis of plotted CC density matrix (t1 to MO)
- Fixed batching estimate in Jacobian D1 doubles
- Moved opening of time-dependent CC files out of constructor
- Fixed copying of HF energy when initializing CC wavefunction
- Removed memory leak in visualization class
- Fixed diagonal scaling factor in CC3 (affects transition moments)
- Sanity checks in Davidson CC ES solver, and changed default for linear dependency threshold that could cause issues in rare cases
- Cluster amplitudes now stored to disk in each iteration of the Newton-Raphson ground state solver
- Fixed Gram-Schmidt biorthonormalization bug that may appear in rare cases for degenerate states (not observed in practice)
- Fixed memory leaks in visualization\_tool, some of which were intel-specific
- Fixed print for oscillator strenths, they are now unitless
- Fixed uninitialized variable field%separation when not specified in input
- Added missing frozen-core contributions to CC dipoles and quadrupoles

### Features
- Integral screening in MLHF for inactive density
- MLCCSD GS, right ES with NTO, CNTO, PAO and Cholesky orbitals and CVS
- Non-linear Davidson solver for CC2 and CC3
- Now possible to do only Cholesky decomposition
- Visualization of CC ground state density and transition densities
- Restart functionality was added to MLCCSD.
- Convention enforced that largest element of an eigenvector from LAPACK has a positive sign.
- Estimate for the maximum amount of memory used printed at the end of the execution.
- Added CMake command to CMake output
- Now possible to restart singles and doubles from singles, and vice versa (e.g. CCS from CCSD and vice versa)
- Now possible to do global restart instead of specifying restart for each individual solver
- Added restart left from right and vice vers with the correct transformation of the basis
- Now possible to run calculations beyond Calcium
- Changes to default thresholds to avoid overly tight thresholds unless requested
- Improved help messages from setup script and added documentation for fortran flags
- Sign of converged Hartree-Fock orbitals are deterministic.
- Minimal Python version is 3.6
- Added additional authors and a reference to the paper to the output
- Added compilation information to the beginning of the output
- Restart added for CNTO orbitals

### Optimization
- Removed unnecessary reordering in Jacobian doubles
- Batching functionality was added to the construction of CNTOs in multilevel coupled cluster.
- Added batching to CCS energy calculation to reduce memory requirement
- Introduced batching in DIIS and Davidson and record storers to minimize copy-operations
- Cholesky vectors are now stored according to occupied-virtual blocks to improve copy/write/read performance
- Reduced prefactor of CC2 and MLCC2 ground state and Jacobian transformation equations by use of the Cholesky vector expression rather than ERIs.
- Reduced prefactor of CC3 Jacobian transpose transformation equations by use of the Cholesky vector expression.
- Reduced prefactor of CCSD and MLCCSD ground state equations by use of the Cholesky vector expression rather than ERIs.
- Removed unnecessary reorderings from CCSD Jacobian transformation

### Structure
- Wavefunction constructors and init routines
- Implemented stream and direct stream files
- Input file read only once
- Removed direct access files
- Geometry now passed to libint in memory
- Record storers for records in memory or on file implemented and used in Davidson
- Improved error messages in memory manager
- The ground state engine now gives an error message if the selected algorithm is not recognized.
- Cleanup and generalization of the biorthogonalization routine
- The vector csiX in the eom/response modules was renamed to xiX according to the greek letter.
- Some improvements in how the autogenerate interfaces Python script recognizes submodules and deals with ampersand, and some restructuring of the script
- Removed direct\_file object
- Excited states are stored in separate stream files which are stored as an array of files in eT.
- Stream files are used for amplitudes and multipliers.
- Separated out convergence tests for solvers into convergence\_tool
- runtest now has a separate file that defines filters
- Python code follows the black code style
- Additional davidson timers
