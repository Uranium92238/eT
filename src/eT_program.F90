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
  !use scf_diis_solver_class
  !use arh_hf_solver_class
!
  !use eri_cd_solver_class
!
  implicit none
!
!   Method allocatable objects
!
    type(hf), allocatable, target    :: hf_wf
    type(mlhf), allocatable, target  :: mlhf_wf 
!
!   Wavefunction pointer
!
    class(wavefunction), pointer :: wf => null()
!
    integer(i15) :: n_methods
!
!  ::::::::::::::::::::::::::::::::::::::::::::::
!  -::-         Print program banner         -::-
!  ::::::::::::::::::::::::::::::::::::::::::::::
!
!    Prepare input and output file
!
    call input%init('eT.inp', 'sequential', 'formatted')
!
    call output%init('eT.out', 'sequential', 'formatted')
    call disk%open_file(output, 'write', 'rewind')
!
    write(output%unit,'(///t16,a)')    'eT - a coupled cluster program'
    write(output%unit,'(t12,a//)') 'S. D. Folkestad, E. F. Kjønstad, 2017-2018'
    flush(output%unit)
!
    n_methods = get_n_methods()
!
    if (n_methods == 0) then
!
!     Do cholesky?
!
    else
!
      if (requested_method('mlhf')) then
!
        allocate(mlhf_wf)
        wf => mlhf_wf
!
      else
!
        allocate(hf_wf)
        wf => hf_wf
!
      endif
!
    endif
!
    stop
!
   !type(scf_diis_solver) :: roothan_hall_hf_solver
   !type(arh_hf_solver)   :: density_minimization_hf_solver
!
   !type(eri_cd_solver) :: chol_solver
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
   !call engine%solve(wf)
    !call db_engine%solve(wf)
!
  ! call roothan_hall_solver%run(wf)
  ! call density_minimization_hf_solver%run(wf)

  !call roothan_hall_hf_solver%run(wf)
!
!  Finalize the wavefunction
!
 !  call wf%finalize()
!
!  Finalize the Libint integral library
!
!  call finalize_libint()
!
!   call disk%close_file(output)
!
end program eT_program
