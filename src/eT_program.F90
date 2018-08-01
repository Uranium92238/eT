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
!  Initialize memory and disk here
!
   call output%init('eT.out', 'sequential', 'formatted')
   call disk%open_file(output, 'write', 'rewind')
!
   write(*,*)'init libint'
   call initialize_libint()
!
! 	Create an SCF engine and ask it to solve the HF wavefunction
!
   write(*,*)'init hf'
   call wf%initialize()
   write(*,*)'init engine'
   call engine%initialize(wf)
!
  !  call wf%integrals%cholesky_decompose(wf%molecule)
  write(*,*)'solve hf'
   call engine%solve(wf) ! solves Hartree Fock
!
   call wf%finalize()
   call finalize_libint()
!
   call disk%close_file(output)
!
end program eT_program
