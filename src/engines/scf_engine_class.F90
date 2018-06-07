module scf_engine_class
!
!!
!!		SCF engine class module
!!		Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
	use kinds
	use file_class
	use hf_class
	use disk_manager_class
!
	implicit none
!
	type :: scf_engine
!
		real(dp) :: energy_threshold   = 1.0D-6
		real(dp) :: residual_threshold = 1.0D-6
!
		integer(i15) :: max_iterations = 100
!
		logical :: restart
!
	contains
!
		procedure :: initialize => initialize_scf_engine
		procedure :: solve 		=> solve_scf_engine
		procedure :: finalize 	=> finalize_scf_engine
!
	end type scf_engine
!
!
contains
!
!
	subroutine initialize_scf_engine(engine, wf)
!!
!! 	Initialize SCF engine
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
		implicit none
!
		class(scf_engine) :: engine
!
		class(hf) :: wf
!
! 		A dummy routine for now, but this routine should read from
! 		file possible changes in the engines variables (thresholds, restart, etc.)
!
	end subroutine initialize_scf_engine
!
!
	subroutine solve_scf_engine(engine, wf)
!!
!!	 	Solve
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!! 	Routine that solves the reference state equations associated with a
!! 	wavefunction
!!
		implicit none
!
		class(scf_engine) :: engine
!
		class(hf) :: wf
!
! 		Initialize wavefunction (this includes reading the molecular
! 		geometry and basis set, fixing parameters, like the number of
! 		atomic orbitals, and initializing arrays to keep in memory)
!
		call wf%initialize()
!
! 		Initialize engine (read thresholds, restart, etc., from file,
! 		but also ask the wavefunction for the number of parameters to solve
! 		for and related information)
!
		call engine%initialize(wf)
!
! 		Solve the equations!
!
! ......
!
! 		Finalize wavefunction (save converged solution to file, for easy
! 		restart, and deallocate arrays that have been kept in memory)
!
		call wf%finalize()
!
! 		Initialize engine (make final deallocations, and other stuff)
!
		call engine%finalize()
!
	end subroutine solve_scf_engine
!
!
	subroutine finalize_scf_engine(engine)
!!
!! 	Finalize SCF engine
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
		implicit none
!
		class(scf_engine) :: engine
!
! 		A dummy routine for now
!
	end subroutine finalize_scf_engine
!
!
end module scf_engine_class
