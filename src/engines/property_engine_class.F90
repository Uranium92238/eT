module property_engine_class 
!!
!!    Coupled cluster property engine class module 
!!    Written by Josefine H. Andersen, February 2019
!!
   use abstract_engine_class
   use ccs_class
   use eri_cd_class
   use davidson_cc_es_class
   use davidson_cc_ip_class
   use davidson_cvs_cc_es_class
   use davidson_cc_multipliers_class
   use davidson_cc_response_class
   use diis_cc_gs_class
   use diis_cc_es_class
   use diis_cc_multipliers_class
   use cc_property_class
!
   type, extends(abstract_engine) :: property_engine
!
      character(len=100) :: algorithm
      character(len=100) :: es_type
!
   contains
!
      procedure :: prepare                   => prepare_property_engine
      procedure :: run                       => run_property_engine
      procedure :: cleanup                   => cleanup_property_engine
!
      procedure :: determine_es_type         => determine_es_type_property_engine
      procedure :: read_algorithm            => read_algorithm_property_engine
!
   end type property_engine
!
contains
! 
   subroutine prepare_property_engine(engine)
!!
!!    Prepare 
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(property_engine) :: engine
!
      engine%name_       = 'Property engine'
!
!     Set standards and then read if nonstandard
!
      engine%algorithm = 'davidson'
      call engine%read_algorithm()
!
      engine%es_type = 'valence'
      call engine%determine_es_type()
!
   end subroutine prepare_property_engine
!
!
   subroutine run_property_engine(engine, wf)
!!
!!    Run  
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(property_engine) :: engine
!
      class(ccs) :: wf
!
      type(eri_cd), allocatable                    :: eri_chol_solver
      type(diis_cc_gs), allocatable                :: cc_gs_solver
      type(diis_cc_es), allocatable                :: cc_es_solver_diis
! 
      type(davidson_cc_es), allocatable            :: cc_valence_es_solver
      type(davidson_cvs_cc_es), allocatable        :: cc_core_es_solver
      type(davidson_cc_multipliers), allocatable   :: cc_multipliers_davidson
      type(diis_cc_multipliers), allocatable       :: cc_multipliers_diis
      type(davidson_cc_response), allocatable      :: cc_response_solver
!
      class(cc_property), pointer :: cc_property_solver
!
      write(output%unit, '(/t3,a,a)') '- Running ', trim(engine%name_)      
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
!     Ground state solver
!
      allocate(cc_gs_solver)
!
      call cc_gs_solver%prepare(wf)
      call cc_gs_solver%run(wf)
      call cc_gs_solver%cleanup(wf)
!
      deallocate(cc_gs_solver)
!
      if (wf%name_ .ne. 'CCS') then
!
         call wf%integrals%write_t1_cholesky(wf%t1)
         call wf%integrals%can_we_keep_g_pqrs()
!
      endif
!
!     Multiplier equation
!
      if (wf%name_ == 'cc2' .or. engine%algorithm == 'diis') then
!
         allocate(cc_multipliers_diis)
!
         call cc_multipliers_diis%prepare(wf)
         call cc_multipliers_diis%run(wf)
         call cc_multipliers_diis%cleanup(wf)
!
         deallocate(cc_multipliers_diis)
!
      elseif (engine%algorithm == 'davidson') then

         allocate(cc_multipliers_davidson)
!
         call cc_multipliers_davidson%prepare(wf)
         call cc_multipliers_davidson%run(wf)
         call cc_multipliers_davidson%cleanup(wf)
!
         deallocate(cc_multipliers_davidson)
!
      endif
!
!     Prepare for excited state
!
      if (engine%algorithm == 'diis') then
!
         allocate(cc_es_solver_diis)
!
         call cc_es_solver_diis%prepare('right')
         call cc_es_solver_diis%run(wf)
         call cc_es_solver_diis%cleanup()
!
         call cc_es_solver_diis%prepare('left')
         call cc_es_solver_diis%run(wf)
         call cc_es_solver_diis%cleanup()
!
         deallocate(cc_es_solver_diis)
!
      elseif (engine%algorithm == 'davidson') then
!
         if (engine%es_type == 'core') then
!
            allocate(cc_core_es_solver)
!
            call cc_core_es_solver%prepare('right')
            call cc_core_es_solver%run(wf)
!
            call cc_core_es_solver%prepare('left')
            call cc_core_es_solver%run(wf)
!
            call cc_core_es_solver%cleanup()
!
            deallocate(cc_core_es_solver)
!
         else ! es_type = valence
!
            allocate(cc_valence_es_solver)
!
            call cc_valence_es_solver%prepare('right')
            call cc_valence_es_solver%run(wf)
!
            call cc_valence_es_solver%prepare('left')
            call cc_valence_es_solver%run(wf)
!
            call cc_valence_es_solver%cleanup()
!
            deallocate(cc_valence_es_solver)
!
         endif
!
      endif
!
! -- TEST
!
      allocate(cc_response_solver)
!
      call cc_response_solver%prepare(wf)
      call cc_response_solver%run(wf)
      !call cc_response_solver%cleanup()
!
      deallocate(cc_response_solver)
!
!     Properties
!
      allocate(cc_property_solver)
!
      call cc_property_solver%prepare(wf)
      call cc_property_solver%run(wf)
      call cc_property_solver%cleanup(wf)
!
      deallocate(cc_property_solver)
!
   end subroutine run_property_engine
!
!
   subroutine cleanup_property_engine(engine)
!!
!!    Cleanup 
!!    Written by Josefine H. Andersen
!!
      implicit none
!
      class(property_engine) :: engine
!
      write(output%unit, '(/t3,a,a)') '- Cleaning up ', trim(engine%name_)
!
   end subroutine cleanup_property_engine
!
!
   subroutine  determine_es_type_property_engine(engine)
!!
!!    Determine excited state type 
!!    Written by JosefineH. Andersen
!!
      implicit none
!
      class(property_engine) :: engine
!
      character(len=100) :: line
!
      integer :: i, n_keywords
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
               !write(output%unit, '(/t3,a,a)') 'Spectra for core excitations not implemented'
               !stop
!
               engine%es_type = 'core'
               return
!
            elseif (line(1:10) == 'valence excitation' ) then
!
               engine%es_type = 'valence'
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
   end subroutine determine_es_type_property_engine
!
!
   subroutine read_algorithm_property_engine(engine)
!!
!!    Read algorithm
!!    Written by Josefine H. Andersen
!!
      implicit none
!
      class(property_engine), intent(inout) :: engine
!
      character(len=100) :: line
!
      integer :: i, n_records
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
   end subroutine read_algorithm_property_engine
!
!
end module property_engine_class
