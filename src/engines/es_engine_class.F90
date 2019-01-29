module es_engine_class
!!
!!    Coupled cluster ground state engine class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
   use abstract_engine_class
   use ccs_class
   use eri_cd_solver_class
   use davidson_cc_es_solver_class
   use davidson_cc_ip_solver_class
   use davidson_cvs_cc_es_solver_class
   use diis_cc_gs_solver_class
   use diis_cc_es_solver_class
   use diis_cc_multipliers_solver_class
!
   type, extends(abstract_engine) :: es_engine 
!
      character(len=100) :: algorithm 
      character(len=100) :: es_type 
!
   contains 
!
      procedure :: prepare                   => prepare_es_engine
      procedure :: run                       => run_es_engine
      procedure :: cleanup                   => cleanup_es_engine
!
      procedure :: determine_es_type         => determine_es_type_es_engine
      procedure :: read_algorithm            => read_algorithm_es_engine
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
      engine%tag       = 'Excited state engine'
!
!     Set standards and then read if nonstandard
!
      engine%algorithm = 'davidson'
      call engine%read_algorithm()
!
      engine%es_type = 'valence'
      call engine%determine_es_type()
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
      type(eri_cd_solver), allocatable              :: eri_chol_solver
      type(diis_cc_gs_solver), allocatable          :: cc_gs_solver
      type(diis_cc_es_solver), allocatable          :: cc_es_solver_diis 
! 
      type(davidson_cc_es_solver), allocatable, target      ::  cc_valence_es
      type(davidson_cvs_cc_es_solver), allocatable, target  ::  cc_core_es
      type(davidson_cc_ip_solver), allocatable, target      ::  cc_valence_ip
!
      class(davidson_cc_es_solver), pointer :: cc_es_solver
!
      write(output%unit, '(/t3,a,a)') '- Running ', trim(engine%tag)
!
!     Cholesky decomposition 
!
      allocate(eri_chol_solver)
!
      call eri_chol_solver%prepare(wf%system)
      call eri_chol_solver%run(wf%system)
!
      call eri_chol_solver%cholesky_vecs_diagonal_test(wf%system)
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
      if (wf%name .ne. 'CCS') then 
!
         call wf%integrals%write_t1_cholesky(wf%t1)
         call wf%integrals%can_we_keep_g_pqrs()
!
      endif
!
!     Prepare for excited state
!
      if (engine%algorithm == 'diis') then 
!
         allocate(cc_es_solver_diis)
!
         call cc_es_solver_diis%prepare()
         call cc_es_solver_diis%run(wf)
         call cc_es_solver_diis%cleanup()
!
         deallocate(cc_es_solver_diis)
!
      elseif (engine%algorithm == 'davidson') then 
!
         if (engine%es_type == 'core') then
!
            allocate(cc_core_es)
            cc_es_solver => cc_core_es
!
            call cc_es_solver%prepare()
            call cc_es_solver%run(wf)
            call cc_es_solver%cleanup()
!
            cc_es_solver => null()
            deallocate(cc_core_es)
!
         elseif(engine%es_type == 'valence ionized') then
!
            allocate(cc_valence_ip)
            cc_es_solver => cc_valence_ip
!
            call cc_es_solver%prepare()
            call cc_es_solver%run(wf)
            call cc_es_solver%cleanup()
!
            cc_es_solver => null()
            deallocate(cc_valence_ip)
!
         elseif(engine%es_type == 'core ionized') then
!
!           Nothing here yet...
!
         else ! es_type = valence
!
            allocate(cc_valence_es)
            cc_es_solver => cc_valence_es
!
            call cc_es_solver%prepare()
            call cc_es_solver%run(wf)
            call cc_es_solver%cleanup()
!
            cc_es_solver => null()
            deallocate(cc_valence_es)
!
         endif
!
      endif
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
!     Nothing here yet...
!
   end subroutine cleanup_es_engine
!
!
   subroutine determine_es_type_es_engine(engine)
!!
!!    Determine excited state type 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(es_engine), intent(inout) :: engine 
!
      character(len=100) :: line
!
      integer(i15) :: i, n_keywords
!
      if (requested_section('cc excited state')) then
!
         call move_to_section('cc excited state', n_keywords)
!
         do i = 1, n_keywords
!
            read(input%unit, '(a100)') line
            line = remove_preceding_blanks(line)
!
            if (line(1:15) == 'core excitation' ) then
!
               engine%es_type = 'core'
               return
!
            elseif (line(1:10) == 'ionization' ) then
!
               engine%es_type = 'valence ionized'
               return
!
            elseif (line(1:15) == 'core ionization' ) then
!
               engine%es_type = 'core ionized'
               return
!
            endif
!
         enddo
!
      endif
!
      engine%es_type = 'valence'
!
   end subroutine determine_es_type_es_engine
!
!
   subroutine read_algorithm_es_engine(engine)
!!
!!    Read algorithm
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none
!
      class(es_engine), intent(inout) :: engine 
!
      character(len=100) :: line
!
      integer(i15) :: i, n_records
!
      call move_to_section('cc excited state', n_records)
!
      do i = 1, n_records
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (line(1:10) == 'algorithm:') then
!
            engine%algorithm = line(11:100)
            engine%algorithm = remove_preceding_blanks(engine%algorithm)
            return
!
         endif
!
      enddo
!
   end subroutine read_algorithm_es_engine
!
!
end module es_engine_class
