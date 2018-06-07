program eT_program
!
!!
!!                        eT - a coupled cluster program
!!         Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!
   use file_class
   use disk_manager_class
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
   call output%init('eT.out', 'sequential', 'formatted')
   call disk%open_file(output, 'write', 'rewind')
!
! 	Create an SCF engine and ask it to solve the HF wavefunction
!
   call wf%initialize()
!
   call engine%solve(wf)
!
   call wf%finalize()
!
   call disk%close_file(output)
!
end program eT_program
