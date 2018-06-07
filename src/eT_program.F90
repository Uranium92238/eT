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
   type(file) :: input, molecule
!
   character(len=40) :: line
   character(len=40) :: current_basis
!
   integer(i15) :: i = 0, n_atoms = 0
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
   call input%init('eT.inp', 'sequential', 'formatted')
   call disk%open_file(input, 'read')
!
   call molecule%init('molecule.inp', 'sequential', 'formatted')
   call disk%open_file(molecule, 'write', 'rewind')
!
   read(input%unit,'(a)') line
   line = remove_preceding_blanks(line)
!
   n_atoms = 0
!
   do while (trim(line) .ne. 'end geometry')
!
      if (line(1:6) == 'basis:' .or.  &
          line(1:6) == 'Basis:' .or.  &
          line(1:6) == 'BASIS:' ) then
!
         current_basis = trim(line(7:40))
         current_basis = remove_preceding_blanks(current_basis)
!
         read(input%unit,'(a)') line
         line = remove_preceding_blanks(line)
!
         do while (trim(line) .ne. 'end geometry'  .and.  &
                    line(1:6) .ne. 'basis:'        .and.  &
                    line(1:6) .ne. 'Basis:'        .and.  &
                    line(1:6) .ne. 'BASIS:')
!
            n_atoms = n_atoms + 1
!
            read(input%unit,'(a)') line
            line = remove_preceding_blanks(line)
!
         enddo
!
         backspace(input%unit)
!
      endif
!
      
      read(input%unit,'(a)') line
      line = remove_preceding_blanks(line)
      write(output%unit, *)trim(line)
!
   enddo
   write(output%unit, *)n_atoms
!
   call disk%close_file(input)
!
   call disk%close_file(output)
!
!
end program eT_program
