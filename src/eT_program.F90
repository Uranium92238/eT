program eT_program
!
!!
!!                        eT - a coupled cluster program
!!         Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!
  use kinds
  use file_class
  use disk_manager_class
  use memory_manager_class
  use libint_initialization
!
  use wavefunction_class
!
  use hf_class
  use uhf_class
  use mlhf_class
!
  use ccs_class
!
  use io_eT_program
!
  use hf_engine_class
  use gs_engine_class
  use abstract_engine_class
!
  use eri_cd_solver_class
!
  implicit none
!
!   Allocatable system
!
    type(molecular_system), allocatable :: system
!
!   Method allocatable objects
!
    type(hf), allocatable, target    :: hf_wf
    type(uhf), allocatable, target   :: uhf_wf
    type(mlhf), allocatable, target  :: mlhf_wf 
!
    type(ccs), allocatable, target :: ccs_wf
!
!   Wavefunction pointer
!
    class(hf), pointer  :: ref_wf => null()
    class(ccs), pointer :: cc_wf  => null()
!
!   Cholesky decomposition solver 
!
    type(eri_cd_solver), allocatable :: chol_solver
!
!   Engines
!
    type(hf_engine), allocatable :: gs_hf_engine
!
    type(gs_engine), allocatable, target :: gs_cc_engine
!
    class(abstract_engine), pointer :: engine => null()
!
    integer(i15) :: n_methods, i
!
    character(len=40), dimension(:), allocatable :: cc_methods
    character(len=40):: cc_engine  
!
!  ::::::::::::::::::::::::::::::::::::::::::::::
!  -::-         Print program banner         -::-
!  ::::::::::::::::::::::::::::::::::::::::::::::
!
!    Prepare input and output file
!
    call output%init('eT.out', 'sequential', 'formatted')
    call disk%open_file(output, 'write', 'rewind')
!
    call input%init('eT.inp', 'sequential', 'formatted')
    call disk%open_file(input, 'read')
!
    write(output%unit,'(///t16,a)')    'eT - a coupled cluster program'
    write(output%unit,'(t12,a//)') 'S. D. Folkestad, E. F. Kjønstad, 2017-2018'
    flush(output%unit)
!
    call initialize_libint()
!
    n_methods = get_n_methods()
!
!   ::  Hartree-Fock calculation (temporarily also cholesky decomposition) 
!
    if (n_methods == 0) then
!
     if (requested_task('cholesky eri')) then
!
        allocate(system)
        call system%prepare()
!
        allocate(chol_solver)
!
        call chol_solver%initialize(system)
        call chol_solver%run(system)
        call chol_solver%finalize()
!
     else
!
       call output%error_msg('no calculation requested.')
!
     endif
!
    else
!
      if (requested_method('mlhf')) then
!
        allocate(mlhf_wf)
        ref_wf => mlhf_wf
!
        call ref_wf%prepare()
!
      else
!
        if (requested_method('uhf')) then 
!
          allocate(uhf_wf)
          ref_wf => uhf_wf
!
        else ! Assume standard RHF 
!
          allocate(hf_wf)
          ref_wf => hf_wf
!
        endif
!
        call ref_wf%prepare()
!
        allocate(gs_hf_engine)
!
        call gs_hf_engine%prepare()     
        call gs_hf_engine%run(ref_wf)     
        call gs_hf_engine%cleanup()     
!
        deallocate(gs_hf_engine)   
!
      endif
!
    endif
!
!   :: Coupled cluster calculation
!
    if (requested_method('mlhf')) n_methods = n_methods - 1
    if (requested_method('hf'))   n_methods = n_methods - 1
    if (requested_method('uhf'))  n_methods = n_methods - 1
!
    if (n_methods .gt. 0) then
!
      allocate(cc_methods(n_methods))
! 
      call read_cc_methods(n_methods, cc_methods)
      call select_engine(cc_engine)
!
      do i = 1, n_methods
!
!       Determine type of CC method
!
        if (cc_methods(i) == 'ccs') then
!
          allocate(ccs_wf)
          cc_wf => ccs_wf
!
        elseif (cc_methods(i) == 'mp2') then
!
        elseif (cc_methods(i) == 'cc2') then
!
        elseif (cc_methods(i) == 'ccsd') then
!
        endif
!
!       Determine engine
!
        if (cc_engine == 'ground state') then
!
          write(output%unit, *) 'hellaw'
!
          allocate(gs_cc_engine)
          engine => gs_cc_engine  
!
        elseif (cc_engine == 'excited state') then
!
        endif
!
!       Solve cc problem
!
        call cc_wf%prepare(ref_wf)
        call ref_wf%cleanup()
        deallocate(ref_wf)
!
        call engine%prepare()
        call engine%run(cc_wf)
        call engine%cleanup()
!
        deallocate(engine)
!
        call cc_wf%cleanup()
        deallocate(cc_wf)
!
      enddo
!
    endif
!
    call finalize_libint()
!
    write(output%unit, '(/t3,a)') 'eT terminated successfully!'
!
   call disk%close_file(output)
   call disk%close_file(input)
!
end program eT_program
