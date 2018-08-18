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
!
  use wavefunction_class
  use hf_class
  use mlhf_class
!
  use io_eT_program
!
  use scf_diis_solver_class
  use arh_hf_solver_class
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
!   Cholesky decomposition solver 
!
    type(eri_cd_solver), allocatable :: chol_solver
!
!   Wavefunction pointer
!
    class(hf), pointer :: wf => null()
!
!   Solvers 
!
    type(scf_diis_solver) :: roothan_hall_hf_solver
    type(arh_hf_solver)   :: density_minimization_hf_solver
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
    n_methods = get_n_methods()
!
    if (n_methods == 0) then
!
     if (requested_task('cholesky eri')) then
!
        allocate(system)
        call system%initialize()
!
        allocate(chol_solver)
!
        call initialize_libint()
!
        call chol_solver%initialize(system)
        call chol_solver%solve(system)
        call chol_solver%finalize()
!
        call finalize_libint()
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
        wf => mlhf_wf
!
        call initialize_libint()
!
        call wf%initialize()
        call wf%finalize()
!
        call finalize_libint()
!
      else
!
        allocate(hf_wf)
        wf => hf_wf
!
        call initialize_libint()
!
        call wf%initialize()
!
        !call roothan_hall_solver%run(wf)
        call density_minimization_hf_solver%run(wf)
        
!
        call wf%finalize()
!
        call finalize_libint()
!
      endif
!
    endif
!
   call disk%close_file(output)
   call disk%close_file(input)
!
end program eT_program

!
   !
   !
!
   !
!
!  Initialize Libint integral library
!
   !call initialize_libint()
!
!  Initialize wavefunction
!
   !call wf%initialize()

  ! call wf%eri_decomp_test_w_active_dens()
!
  !  call chol_solver%initialize(wf%system)
  !  call chol_solver%solve(wf%system)
  !  call chol_solver%finalize()
!
!  Ask the Hartree-Fock (HF) solver to find the HF solution
!
!  Finalize the wavefunction
!
 !  call wf%finalize()
!
!  Finalize the Libint integral library
!
!  call finalize_libint()
