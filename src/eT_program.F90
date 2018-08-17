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
   use io_utilities
!
   use hf_class
   use mlhf_class
!
   use scf_diis_solver_class
   use arh_hf_solver_class
!
   use eri_cd_solver_class
!
   implicit none
!
   type(scf_diis_solver) :: roothan_hall_hf_solver
   type(arh_hf_solver)   :: density_minimization_hf_solver
!
   type(eri_cd_solver) :: chol_solver
!
   type(hf) :: wf
  ! type(mlhf) :: wf
!
!
!  Initialize memory and disk here
!
   call output%init('eT.out', 'sequential', 'formatted')
   call disk%open_file(output, 'write', 'rewind')
!
!  ::::::::::::::::::::::::::::::::::::::::::::::
!  -::-         Print program banner         -::-
!  ::::::::::::::::::::::::::::::::::::::::::::::
!
   write(output%unit,'(///t16,a)')    'eT - a coupled cluster program'
   write(output%unit,'(t12,a//)') 'S. D. Folkestad, E. F. Kjønstad, 2017-2018'
   flush(output%unit)
!
!  Initialize Libint integral library
!
   call initialize_libint()
!
!  Initialize wavefunction
!
   call wf%initialize()
  ! call wf%eri_decomp_test_w_active_dens()
!
  !  call chol_solver%initialize(wf%system)
  !  call chol_solver%solve(wf%system)
  !  call chol_solver%finalize()
!
!  Ask the Hartree-Fock (HF) solver to find the HF solution
!
   call density_minimization_hf_solver%run(wf)
!
!  Finalize the wavefunction
!
   call wf%finalize()
!
!  Finalize the Libint integral library
!
   call finalize_libint()
!
   call disk%close_file(output)
!
end program eT_program
