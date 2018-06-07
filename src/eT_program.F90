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
   use hf_engine_class
!
   implicit none
!
   type(hf_engine) :: engine
!
   type(hf) :: wf
!
   type(file) :: input
!
   character(len=40) :: line
!
   integer(i15) :: i = 0
!
!  Initialize memory and disk here
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
  ! call input%init('eT.inp', 'sequential', 'formatted')
  ! call disk%open_file(input, 'read')
  ! do i = 1, 5
  !    read(input%unit) line
  !    write(output%unit,*) check_if_blank(line)
  ! enddo
  ! call disk%close_file(input)
!
!
end program eT_program
