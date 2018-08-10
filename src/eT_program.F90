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
   use integral_manager_class
!
   use hf_class
   use mlhf_class
   use hf_engine_class
   use dmm_hf_engine_class
!
   use eri_cd_engine_class
!
   implicit none
!
   type(hf_engine) :: engine

   type(dmm_hf_engine) :: density_hf_engine
   type(eri_cd_engine) :: chol_engine
!
   type(hf) :: wf
   !type(mlhf) :: wf
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
  !  call chol_engine%initialize(wf%system)
  !  call chol_engine%solve(wf%system)
  !  call chol_engine%finalize()
!!
!  Ask the Hartree-Fock (HF) engine to find the HF solution
!
 !  call engine%solve(wf)
   call density_hf_engine%solve(wf)
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
