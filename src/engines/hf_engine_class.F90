module hf_engine_class
!!
!!    Hartree-Fock engine class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
   use hf_class 
!
   use dmm_hf_solver_class
   use hf_solver_class
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
      type(dmm_hf_solver), allocatable, target :: density_minimization_hf_solver
      type(hf_solver), allocatable, target     :: roothan_hall_hf_solver
!
      class(dmm_hf_solver), pointer :: solver => null() ! Change class to common ancestor when this is done 
!
      allocate(density_minimization_hf_solver)
      solver => density_minimization_hf_solver
!
      call solver%run(wf)
!
      deallocate(density_minimization_hf_solver)
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
