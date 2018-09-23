module es_engine_class
!!
!!    Coupled cluster ground state engine class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
   use abstract_engine_class
   use ccs_class
   use eri_cd_solver_class
   use davidson_cc_es_solver_class
   use diis_cc_gs_solver_class
!
   type, extends(abstract_engine) :: es_engine 
!
   contains 
!
      procedure :: prepare => prepare_es_engine
      procedure :: run     => run_es_engine
      procedure :: cleanup => cleanup_es_engine
!
      procedure :: print_banner    => print_banner_es_engine
      procedure :: print_summary   => print_summary_es_engine
!
   end type es_engine 
!
contains
!
   subroutine prepare_es_engine(engine)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(es_engine) :: engine 
!
   end subroutine prepare_es_engine
!
!
   subroutine run_es_engine(engine, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(es_engine) :: engine 
!
      class(ccs) :: wf
!
      type(eri_cd_solver), allocatable          :: eri_chol_solver
      type(diis_cc_gs_solver), allocatable      :: cc_gs_solver 
      type(davidson_cc_es_solver), allocatable  :: cc_es_solver 
!
!     Cholesky decoposition 
!
      allocate(eri_chol_solver)
!
      call eri_chol_solver%prepare(wf%system)
      call eri_chol_solver%run(wf%system)
!
      call eri_chol_solver%cholesky_vecs_diagonal_test(wf%system)
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
!     Prepare for excited state
!
      call wf%initialize_t1()
      wf%t1 = zero
!
!     Excited state solution 
!
      allocate(cc_es_solver)
!
      call cc_es_solver%prepare(wf)
      call cc_es_solver%run(wf)
      call cc_es_solver%cleanup(wf)
!
      deallocate(cc_es_solver)
!
   end subroutine run_es_engine
!
!
   subroutine cleanup_es_engine(engine)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(es_engine) :: engine 
!
   end subroutine cleanup_es_engine
!
!
   subroutine print_banner_es_engine(engine)
!!
!!    Print banner 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(es_engine) :: engine 
!
   end subroutine print_banner_es_engine
!
!
   subroutine print_summary_es_engine(engine)
!!
!!    Print summary
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(es_engine) :: engine 
!
   end subroutine print_summary_es_engine
!
!
end module es_engine_class
