# eT v1.9.0
### Features
- Added total and iteration timers to solvers. eT-program/eT!1109
- Added option to zero out arrays when allocating. eT-program/eT!1112
- Added option to set logical arrays to true/false when allocating. eT-program/eT!1128

### Structure
- Added separate engine for Lanczos and removed `cc_es_eigen_davidson_print_tool`. eT-program/eT!1095
- Removed unused routines in CCS and MLCCSD. eT-program/eT!1105
- Made all unformatted files stream files and removed unformatted sequential. eT-program/eT!1118, eT-program/eT!1132
- Restructured I/O classes moving more functionality into `abstract_file_class`. eT-program/eT!1118, eT-program/eT!1132
- Removed unused `abstract_hf_solver_class.F90`. eT-program/eT!1121
- Removed warnings from intel compilers. eT-program/eT!1128
- Response equations in CC theory are now solved using general Davidson solver. eT-program/eT!1133

### Tests
- Added checks for cube and auxiliary files to the time-dependent CC tests. eT-program/eT!1096
- Removed for loop over the inputs in the test scripts. eT-program/eT!1098
- Frozen HF tests with excitation energies and transition moments. eT-program/eT!1101
- Added unit test for z-matrix tool. eT-program/eT!1103
- Added restart test for CCS and 3-level MLCCSD tests. eT-program/eT!1105
- Simplified the use of keyword arguments in `filters.py`. eT-program/eT!1106

### CI
- Updated version of `fortran-code-quality` to v1.3.0. eT-program/eT!1139
- CI now always runs python and fortran code-quality, but stops if one of these stages does not succeed. eT-program/eT!1140

# eT v1.8.8
### Bugfixes
- Fixed bug where libint threw an error because of whitespace in the final lines of a basis set file by removing trailing whitespace from g94 files. eT-program/eT!1131

>>>>>>> 6d5a14c30c9ce9e02657ccf7f2ecbd8bc14bb56d
# eT v1.8.7
### Bugfixes
- Revert update of Libint library, as the eT compilation fails. eT-program/eT!1125
- Fixed debug print in `memory_tracker.F90` that crashed because it expected integers but received strings (eT-program/eT#586). eT-program/eT!1125
- Fixed segfault that appears in EOM properties with Intel 20.2.6 (eT-program/eT#586). eT-program/eT!1125

# eT v1.8.6
### Features
- Update of Libint library that ships with eT to include higher angular momentum. eT-program/eT!1119

# eT v1.8.5
### Bugfixes
- Fixed bug in the FCI linear transformation where $`E_{qp}`$ was applied instead of $`E_{pq}`$ . eT-program/eT!1115

# eT v1.8.4
### Bugfixes
- Fixed bug where visualizing the active hf density was giving a plot of the hf density. eT-program/eT!1113

# eT v1.8.3
### Bugfixes
- Fixed uninitialized `exists_frozen_fock_terms` in FCI. eT-program/eT!1110

# eT v1.8.2
### Bugfixes
- Fix construction of virtual Cholesky orbitals in 3-level MLCCSD. eT-program/eT!1107

# eT v1.8.1
### Bugfixes
- Fixes bug where keyword `coulomb exchange terms` was not read in MLHF. eT-program/eT!1104

# eT v1.8.0
### Features
- Added calculation of FCI densities and dipole and quadrupole tests. eT-program/eT!1087
- Frequency-dependent (and static) polarizabilities for HF. eT-program/eT!1054

### Structure
- Changed routines in parameters.F90 to increase covarage. eT-program/eT!1068
- Removed unused routines in stream\_file.F90. eT-program/eT!1069
- Removed unused routines in array\_utilities.F90. eT-program/eT!1071
- Restructured eT program and added factories for engines and wavefunctions. eT-program/eT!1074 and eT-program/eT!1080
- Remove frequency from transformation tool. Dummy frequencies removed from general linear and general eigen davidson solvers. eT-program/eT!1076
- Restructured eT program and added factories for engines and wavefunctions. eT-program/eT!1074
- `Transformation_tools` renamed to transformation. eT-program/eT!1054
- Davidson tool now counts the number of new trials. eT-program/eT!1089
- Introduced tools to build the G-matrix (two-electron contribution to Fock) or Coulomb and exchange matrices separately. eT-program/eT!870
- Restructured printing of mean values. eT-program/eT!1086
- Restructured FCI engine to use tasks. eT-program/eT!1087
- Restructured reference/HF engines to use task structure. eT-program/eT!1088
- Added factories for cc start vector and cc projection tools. eT-program/eT!1091
- Added storers that know how to store amplitudes. eT-program/eT!1092 eT-program/eT!1094

### Optimization
- Removed N^3 scaling from Fock construction visible for very large systems. eT-program/eT!870
- Linear exchange algorithm added (similar to LinK). eT-program/eT!870
- Non-scaling iterative Fock matrix construction and linear scaling non-iterative Fock matrix construction for MLHF. eT-program/eT!870

# eT v1.7.4
### Bugfixes
- Fixes problem with the use of bath orbitals when cc-in-hf (frozen hf) is used. Adding test with valence ionization and frozen hf. eT-program/eT!1084

# eT v1.7.3
### Bugfixes
- Fixed bug where python version 3.10 was wrongly recognized to be older than 3.6. eT-program/eT!1079

# eT v1.7.2
### Bugfixes
- Fixed bug in spin-multiplicities, where sum of squares was incorrectly assumed equal to square of sum. Happened only with -O0 flag. eT-program/eT!1072


# eT v1.7.1
### Bugfixes
- Fixed OpenMP bug where thread-integer was not private. Happened only with -O0 flag. eT-program/eT!1066

# eT v1.7.0
### Bugfixes
- QED keyword photons removed. eT-program/eT!1051

### Features
- Linear response now available at the CC2 level (polarizabilities and transition strengths). eT-program/eT!1039
- Restart for the Lanczos solver. eT-program/eT!1036
- FCI for ground and excited states of closed and open-shell systems. eT-program/eT!1060

### Tests
- Added unit tests for for untested routines in `tools/maps`. eT-program/eT!1035
- Added for plotting of the active hf density. eT-program/eT!1045

### Structure
- Cleanup of engines, by introducing tasks that are called from the engine. eT-program/eT!968, eT-program/eT!1040, eT-program/eT!1044
- Using the eigen_davidson_solver for the CC excited state equations instead of the cc_es_davidson_solver. The cc_es_davidson_solver is deleted. eT-program/eT!1038, eT-program/eT!1050
- Projection tools now ask CC wave functions to project out unwanted components. The projection vectors have been deleted. eT-program/eT!1048
- eT stops with error if there are no electrons. eT-program/eT!1061

# eT v1.6.4
### CI
- Updated black to version 22.3.0 because a dependency of black was updated making black 22.1.0 crash. eT-program/eT!1057

# eT v1.6.3
### Bugfixes
- Fixed integer overflow in `batch_setup`. The product of two `int32`
overflowed before being converted into `int64`. eT-program/eT!1053

# eT v1.6.2
### Bugfixes
- Fixed visualization with cartesian basis sets. eT-program/eT!1043

# eT v1.6.1
### README
- Updated README with information for libint and pcmsolver. eT-program/eT!1041

# eT v1.6.0
### Features
- Print tool for eigen and linear Davidson solvers. eT-program/eT!983
- Added timers for calculate energy. eT-program/eT!993
- Memory tracker now prints the batching tag and memory difference. eT-program/eT!1000
- References are now collected in eT_references.F90. eT-program/eT!1019
- Added timers for visualization. eT-program/eT!1024

### Optimization
- Coupled cluster Fock matrix construction is N^4 scaling. eT-program/eT!973
- Reduce memory usage for MLCC energy calculation (reusing CCS routine). eT-program/eT!974
- Optimized omega lowmemory CC2 using a similar algorithm as for CC3. eT-program/eT!986
- New Cholesky tools can now load blocks into memory to be used in e.g. integral construction. eT-program/eT!1015

### Structure
- New ERI and Cholesky vector tools. Autogeneration of complex tools. eT-program/eT!902
- Cleanup of use statements in solvers, molecule and observer. eT-program/eT!1012
- Fix typo in `construct_mo_basis_transformation`. eT-program/eT!1020
- Removed unused routines and added unit tests for `range_class` and `memory_manager`. eT-program/eT!1022
- Added abstract solver class. eT-program/eT!1018

### Tests
- Code quality pipeline now detects unexpected indentation level in Fortran source files. eT-program/eT!1017
- QED-HF tests are now run with default print level. eT-program/eT!1023

# eT v1.5.15
### Bugfixes
- Fixed bug where LR polarizabilities would be incorrect in cases where non-diagonal contributions are requested without the associated diagonal contributions. eT-program/eT!1025
### CI
- Now using Black version 22.1.0. eT-program/eT!1028

# eT v1.5.14
### Tests
- Intel CI stage (like GNU) split into two jobs: build and run. eT-program/eT!1014

# eT v1.5.13
### Tests
- Intel CI stage now runs in merge requests, and compressed the build and release stages. eT-program/eT!1009

# eT v1.5.12
### Tests
- Added stage to test Intel compilers (ifort/icpc) in Gitlab CI/CD. eT-program/eT!1004

# eT v1.5.11
### Bugfixes
- Fixed format in call to printf in `eri_cd_class.F90`. eT-program/eT!1001

# eT v1.5.10
### Bugfixes
- Fixed memory estimate in `construct_cholesky_t1_oo`. eT-program/eT!998

# eT v1.5.9
### Tests
- Fix gcov call for MKL-Ubuntu Docker image. eT-program/eT!996

# eT v1.5.8
### Tests
- Fix verification of codecov uploader for MKL-Ubuntu Docker image. eT-program/eT!995

# eT v1.5.7
### Tests
- Updated CI/CD setup to use MKL-Ubuntu Docker image and adapt to dedicated eT test machine. eT-program/eT!992

# eT v1.5.6
### License
- Updated Readme and added license to `index_invert`. eT-program/eT!990

# eT v1.5.5
### Copyright
- Updated copyright year to 2022. eT-program/eT!988

# eT v1.5.4
### Bugfixes
- Circumvented intel segmentation violation by removing use statements for types in submodules again. eT-program/eT!987

# eT v1.5.3
### Bugfixes
- Fixed bugs that appear at -O0 optimization level, and added -O0 to coverage pipelines. eT-program/eT!982

# eT v1.5.2
### Bugfixes
- Coverage files produced for codecov code coverage. eT-program/eT!979

# eT v1.5.1
### Bugfixes
- Update submodules again. eT-program/eT!977

# eT v1.5.0

### Bugfixes
- Two workarounds added to circumvent bugs in the Intel compiler (ifort). File destructors no longer check whether the file is open. Two OpenMP loops are incorrectly optimized with -O3 and have been modified to circumvent incorrect compiler optimization. eT-program/eT!936
- Disk storage default for davidson solver to obtain multipliers. eT-program/eT!942
- `max_dim_red` of excited state Davidson solvers is now set to `max(100,10*n_singlet_states)` after the number of singlet states have been read from input. eT-program/eT!932

### Features
- Using orthogonal AO (OAO) basis for HF/UHF gradient in Roothan-Hall SCF solver. eT-program/eT!914
- When the code is batching a tag will be printed now for verbose print level. eT-program/eT!919
- Now possible to specify start guesses for excited states as `state guesses: {i=1,a=2}, {i=1,a=3}, ...`. eT-program/eT!933
- Convergence testing now with tolerance of +-1, to avoid fails due to numerics. eT-program/eT!938
- Added CI pipeline that runs a simple program to parse the Fortran code and identify various code-issues. eT-program/eT!937
- Added check for close lying/overlapping atoms. eT-program/eT!943
- Added new reference method: QED Hartree-Fock (qed-hf). eT-program/eT!941
- Added basis sets ccX-nZ basis sets. eT-program/eT!955
- ROHF and CUHF. eT-program/eT!957
- Cartesian BFGS geometry optimization solver (for Hartree-Fock theory) replaced by a more robust BFGS solver that uses redundant internal coordinates. eT-program/eT!918
- Multimodel CCSD/CC2 is now available for the ground state amplitudes. eT-program/eT!950

### Structure
- Moved the SCF preparations for HF outside of SCF solver. eT-program/eT!908
- Merged the `abstract_convergence_tool` and the `convergence_tool` classes. eT-program/eT!909
- SAD handled in the sad\_tool\_class. eT-program/eT!926
- Removed unused routines in reordering. eT-program/eT!939
- Moved a Fock matrix term from CCSD Jacobian into doubles Jacobian. eT-program/eT!949
- G(D) definition changed to (2g\_wxyz - g\_wzyx)D\_yz.  eT-program/eT!964
- Cleanup of use statements in HF hierarchy. eT-program/eT!966
- Removed unused routine in abstract\_file\_class.F90. eT-program/eT!967
- Cleanup of use statements in CC hierarchy. eT-program/eT!969
- No module imports all reordering and array utilities anymore. eT-program/eT!970
- Added unit tests for angular momentum class. eT-program/eT!972

### Optimization
- One-electron and effective contributions to the CC-Fock matrix are now calculated once instead of multiple times per Fock construction. eT-program/eT!925
- Removed an unnecessary non-iterative o3v3 term (the G2-2 intermediate) as well as an iterative o3v3 term in CCSD Jacobian transpose G2. eT-program/eT!946
- Removed two o3v3 terms and an unnecessary g\_ovvv construction in CCSD Jacobian transpose G1. eT-program/eT!944
- Removed an o3v3 term in CCSD Jacobian transpose B2. eT-program/eT!945
- Removed an unnecessary o3v3 term in CCSD Jacobian transpose C2. eT-program/eT!946
- Removed an o3v3 term in Jacobian CCSD D2 and reduced integral costs in the term. eT-program/eT!951
- Replaced an o1v3nJ by an o2v2nJ contraction in Jacobian transpose doubles B1. eT-program/eT!948
- Removed an iterative and a non-iterative o3v3 term in CCSD Jacobian H2. eT-program/eT!958
- Pipeline runs every n-th test (n: number of threads) starting with test i (i: thread number). eT-program/eT!960
- Reduced size of Lanczos test for valence states. eT-program/eT!960
- Removed an o3v3 term in CCSD Jacobian I2. eT-program/eT!959

# eT v1.4.4
### Bugfixes
- Find substring in string also if it contains the last character. eT-program/eT!952

# eT v1.4.3
### Bugfixes
- Update list of authors. eT-program/eT!927

# eT v1.4.2
### Bugfixes
- Fixed bug that made it impossible to enable the forced batching expert option. eT-program/eT!916

# eT v1.4.1
### Bugfixes
- Changed warnings about wrong memory estimates to verbose prints. eT-program/eT!911

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
- Now possible to run linear response calculations at the CCSD level. eT-program/eT!890
- Default algorithm for ground state in CC3 changed to multimodel Newton. eT-program/eT!892
- Full and multimodel Newton-Raphson algorithm now available for ground state multipliers. The default algorithm for CC3 is changed from DIIS to multimodel Newton. eT-program/eT!895
- Jacobian transpose transformation for CCS with Cholesky vectors, scaling is N^4. eT-program/eT!897
- SAD can now run with charges on atoms. eT-program/eT!832

### Optimization
- Adding option for overlap screening to CC calculation of one electron integral. eT-program/eT!851
- Record storers now delete their file when finalized is called. eT-program/eT!893
- Optimization of E1 term in the CCSD Jacobian transpose transformation. eT-program/eT!896

### Structure
- Added block class using the range class. Facilitates the handling of multiple index ranges. eT-program/eT!823
- Angular momentum tools handle components of the angular momenta. eT-program/eT!864
- Generalized davidson solver for linear equations. eT-program/eT!901

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
