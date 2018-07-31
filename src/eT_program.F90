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
   use hf_engine_class
!
   implicit none
!
   type(hf_engine) :: engine
!
   type(hf) :: wf
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
!
!  Ask the Hartree-Fock (HF) engine to find the HF solution
!
   call engine%solve(wf)
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
