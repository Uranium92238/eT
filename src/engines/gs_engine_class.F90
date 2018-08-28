module gs_engine_class
!!
!!    Hartree-Fock engine class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
   use abstract_engine_class
   use ccs_class
!
   type, extends(abstract_engine) :: gs_engine 
!
   contains 
!
      procedure :: prepare => prepare_gs_engine
      procedure :: run     => run_gs_engine
      procedure :: cleanup => cleanup_gs_engine
!
      procedure :: print_banner    => print_banner_gs_engine
      procedure :: print_summary   => print_summary_gs_engine
!
!
   end type gs_engine 
!
contains
!
   subroutine prepare_gs_engine(engine)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(gs_engine) :: engine 
!
   end subroutine prepare_gs_engine
!
!
   subroutine run_gs_engine(engine, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(gs_engine) :: engine 
!
      class(ccs) :: wf
!
   end subroutine run_gs_engine
!
!
   subroutine cleanup_gs_engine(engine)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(gs_engine) :: engine 
!
   end subroutine cleanup_gs_engine
!
   subroutine print_banner_gs_engine(engine)
!!
!!    Print banner 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(gs_engine) :: engine 
!
   end subroutine print_banner_gs_engine
!
   subroutine print_summary_gs_engine(engine)
!!
!!    Print summary
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(gs_engine) :: engine 
!
   end subroutine print_summary_gs_engine
!
end module gs_engine_class
