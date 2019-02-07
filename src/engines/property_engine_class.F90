module property_engine_class 
!!
!!    Coupled cluster property engine class module 
!!    Written by Josefine H. Andersen, February 2019
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
! -- here allocatable variables (solvers)
!
      write(output%unit, '(/t3,a,a)') '- Running ', trim(engine%name_)      
!
!     Cholesky decomposition
!

!
!     Ground state solver
!

!
!     Prepare for excited state
!

!
!     Multiplier equation
!

!
   end run_property_engine
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
      if (requested_section('cc property')) then
!
         call move_to_section('cc property', n_keywords)
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
   end subroutine determine_es_type_property_engine(engine)
!
!
   subroutine read_algorithm_es_engine(engine)
!!
!!    Read algorithm
!!    Written by Josefine H. Andersen
!!
      implicit none
!
      class(es_engine), intent(inout) :: engine
!
      character(len=100) :: line
!
      integer :: i, n_records
!
      call move_to_section('cc property', n_records)
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
end module property_engine_class
