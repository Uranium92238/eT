module gs_engine_class
!!
!!    Coupled cluster ground state engine class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
   use abstract_engine_class
   use ccs_class
   use eri_cd_solver_class
   use diis_cc_gs_solver_class
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
      type(eri_cd_solver), allocatable     :: eri_chol_solver
      type(diis_cc_gs_solver), allocatable :: cc_gs_solver 
!
!     Cholesky decoposition 
!
      allocate(eri_chol_solver)
!
      call eri_chol_solver%prepare(wf%system)
      call eri_chol_solver%run(wf%system)
!
      call eri_chol_solver%cholesky_vecs_diagonal_test()
      call eri_chol_solver%full_test_cholesky_vecs(wf%system)
!
      call eri_chol_solver%construct_mo_cholesky_vecs(wf%system, wf%n_mo, wf%orbital_coefficients)
!
      call wf%integrals%prepare(eri_chol_solver%n_cholesky, wf%n_o, wf%n_v)
!
      call eri_chol_solver%cleanup()
      deallocate(eri_chol_solver)
!
!     Ground state solution 
!
      allocate(cc_gs_solver)
!
      call cc_gs_solver%prepare(wf)
      call cc_gs_solver%run(wf)
      call cc_gs_solver%cleanup(wf)
!
      deallocate(cc_gs_solver)
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
!
end module gs_engine_class
