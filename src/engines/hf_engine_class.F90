module hf_engine_class
!!
!!    Hartree-Fock engine class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
   use hf_class 
!
   use arh_hf_solver_class
   use scf_diis_solver_class
   use scf_solver_class
!
   type hf_engine 
!
   contains 
!
      procedure :: prepare    => prepare_hf_engine
      procedure :: run        => run_hf_engine
      procedure :: cleanup    => cleanup_hf_engine
!
      procedure :: read_algorithm => read_algorithm_hf_engine
!
   end type hf_engine 
!
contains
!
!
   subroutine prepare_hf_engine(engine)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(hf_engine) :: engine 
!
   end subroutine prepare_hf_engine
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
      type(arh_hf_solver), allocatable   :: arh
      type(scf_diis_solver), allocatable :: scf_diis
      type(scf_solver), allocatable      :: scf
!
      character(len=100) :: algorithm
!
      if (requested_section('hf')) then
!
         call engine%read_algorithm(algorithm)
!
         if (trim(algorithm) == 'aug-rh') then
!
            call output%error_msg("aug-rh is not yet working (to-do)") 
!
            allocate(arh)
!
            call arh%prepare(wf)
            call arh%run(wf)
            call arh%cleanup(wf)
!
            deallocate(arh)
!
         elseif (trim(algorithm) == 'scf-diis') then
!
            allocate(scf_diis)
!
            call scf_diis%prepare(wf)
            call scf_diis%run(wf)
            call scf_diis%cleanup(wf)
!
            deallocate(scf_diis)
!
         elseif (trim(algorithm) == 'scf') then 
!
            allocate(scf)
!
            call scf%prepare(wf)
            call scf%run(wf)
            call scf%cleanup(wf)
!
            deallocate(scf)
!
         else
!
            call output%error_msg('did not recognize hf algorithm: '// algorithm)
!
         endif
!
      else ! Default: use SCF DIIS algorithm 
!
         allocate(scf_diis)
!
         call scf_diis%prepare(wf)
         call scf_diis%run(wf)
         call scf_diis%cleanup(wf)
!
         deallocate(scf_diis)
!
      endif
!
   end subroutine run_hf_engine
!
!
   subroutine cleanup_hf_engine(engine)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(hf_engine) :: engine 
!
   end subroutine cleanup_hf_engine
!
!
   subroutine read_algorithm_hf_engine(engine, algorithm)
!!
!!    Read algorithm
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none
!
      class(hf_engine), intent(in) :: engine 
!
      character(len=100), intent(out) :: algorithm
!
      character(len=100) :: line
!
      integer(i15) :: i, n_records
!
      call move_to_section('hf', n_records)
!
      do i = 1, n_records
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (line(1:10) == 'algorithm:') then
!
            algorithm = line(11:100)
            algorithm = remove_preceding_blanks(algorithm)
            return
!
         endif
!
      enddo
!
      algorithm = 'scf-diis' ! Standard
!
   end subroutine read_algorithm_hf_engine
!
!
end module hf_engine_class
