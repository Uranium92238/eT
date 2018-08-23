module hf_engine_class
!!
!!    Hartree-Fock engine class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
   use hf_class 
!
   use hf_solver_class
   use arh_hf_solver_class
   use scf_diis_solver_class
!
   type hf_engine 
!
   contains 
!
      procedure :: initialize => initialize_hf_engine
      procedure :: run        => run_hf_engine
      procedure :: finalize   => finalize_hf_engine
!
   end type hf_engine 
!
contains
!
   subroutine initialize_hf_engine(engine)
!!
!!    Initialize 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(hf_engine) :: engine 
!
   end subroutine initialize_hf_engine
!
!
   subroutine run_hf_engine(engine, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(hf_engine) :: engine 
!
      class(hf) :: wf 
!
      type(arh_hf_solver), allocatable, target   :: arh_solver 
      type(scf_diis_solver), allocatable, target :: scf_solver
!
      class(hf_solver), pointer :: solver => null()
!
      ! allocate(arh_solver)
      ! solver => arh_solver
!
      allocate(scf_solver)
      solver => scf_solver
!
      call solver%run(wf)
!
     ! deallocate(arh_solver)
      deallocate(scf_solver)
!
   end subroutine run_hf_engine
!
!
   subroutine finalize_hf_engine(engine)
!!
!!    Finalize 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(hf_engine) :: engine 
!
   end subroutine finalize_hf_engine
!
!
end module hf_engine_class
