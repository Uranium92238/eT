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
      type(arh_hf_solver), allocatable, target   :: arh
      type(scf_diis_solver), allocatable, target :: scf_diis
      type(scf_solver), allocatable, target      :: scf
!
      class(hf_solver), pointer :: solver => null()
!
      character(len=100) :: algorithm
!
      if (requested_section('hf')) then
!
         call engine%read_algorithm(algorithm)
!
         if (trim(algorithm) == 'aug-rh') then
!
            allocate(arh)
            solver => arh
!
            call solver%prepare(wf)
            call solver%run(wf)
            call solver%cleanup(wf)
!
            deallocate(arh)
!
         elseif (trim(algorithm) == 'scf-diis') then
!
            select type (wf)
            type is (hf) 
!
               allocate(scf_diis)
               solver => scf_diis
!
               call solver%prepare(wf)
               call solver%run(wf)
               call solver%cleanup(wf)
!
               deallocate(scf_diis)
!
            class default
!
               call output%error_msg('SCF DIIS does not support wavefunction ' // wf%name)
!
            end select
!
         elseif (trim(algorithm) == 'scf') then 
!
            allocate(scf)
            solver => scf 
!
            call solver%prepare(wf)
            call solver%run(wf)
            call solver%cleanup(wf)
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
         solver => scf_diis
!
         call solver%prepare(wf)
         call solver%run(wf)
         call solver%cleanup(wf)
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
         write(output%unit, *) trim(line)
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
      algorithm = 'scf-diis'
!
   end subroutine read_algorithm_hf_engine
!
!
end module hf_engine_class
