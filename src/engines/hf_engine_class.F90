module hf_engine_class
!
!!
!!		HF engine class module
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
   type :: hf_engine
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
      procedure :: initialize => initialize_hf_engine
      procedure :: solve      => solve_hf_engine
      procedure :: finalize   => finalize_hf_engine
!
   end type hf_engine
!
!
contains
!
!
   subroutine initialize_hf_engine(engine, wf)
!!
!!    Initialize HF engine
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf_engine) :: engine
!
      class(hf) :: wf
!
!     A dummy routine for now, but this routine should read from
!     file possible changes in the engines variables (thresholds, restart, etc.)
!
   end subroutine initialize_hf_engine
!
!
   subroutine solve_hf_engine(engine, wf)
!!
!!    Solve
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Routine that solves the reference state equations associated with a
!!    wavefunction
!!
      implicit none
!
      class(hf_engine) :: engine
!
      class(hf) :: wf
!
!     Initialize engine (read thresholds, restart, etc., from file,
!     but also ask the wavefunction for the number of parameters to solve
!     for and related information)
!
      call engine%initialize(wf)
!
!     Solve the equations!
!
!
!     Initialize engine (make final deallocations, and other stuff)
!
      call engine%finalize()
!
   end subroutine solve_hf_engine
!
!
   subroutine finalize_hf_engine(engine)
!!
!! 	Finalize SCF engine
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf_engine) :: engine
!
!     A dummy routine for now
!
   end subroutine finalize_hf_engine
!
!
end module hf_engine_class
