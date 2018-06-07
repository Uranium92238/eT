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
	use scf_engine_class
!
	use hf_class
!
	implicit none
!
	type(scf_engine) :: engine
!
	type(hf) :: wf
!
	output%name = 'eT.out'
	call disk%open_file(output, 'formatted', 'write', 'sequential')
!
! 	Create an SCF engine and ask it to solve the HF wavefunction
!
	call engine%solve(wf)
!
	call disk%close_file(output)
!
end program eT_program
