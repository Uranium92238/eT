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
  use libint_initialization
!
  use wavefunction_class
!
  use hf_class
  use mlhf_class
!
  use ccs_class
!
  use io_eT_program
!
  use hf_engine_class
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
    type(mlhf), allocatable, target  :: mlhf_wf 
!
    type(ccs), allocatable, target :: ccs_wf
!
!   Wavefunction pointer
!
    class(hf), pointer  :: ref_wf => null()
    class(ccs), pointer :: cc_wf => null()
!
!   Cholesky decomposition solver 
!
    type(eri_cd_solver), allocatable :: chol_solver
!
!   Engines
!
    type(hf_engine), allocatable :: gs_hf_engine
!
    integer(i15) :: n_methods
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
        call chol_solver%solve(system)
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
        allocate(hf_wf)
        ref_wf => hf_wf
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
    if (n_methods .gt. 0) then
!
      if (requested_method('ccs')) then
!
        allocate(ccs_wf)
        cc_wf => ccs_wf
!
        write(output%unit, *)'CCS!'
        flush(output%unit)
!
        call cc_wf%initialize(ref_wf)
        call ref_wf%cleanup()
        deallocate(ref_wf)
!
      endif
!
    endif
!
    call finalize_libint()
!
   call disk%close_file(output)
   call disk%close_file(input)
!
end program eT_program
