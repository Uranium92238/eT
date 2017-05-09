program eT_program
!
!
!                 Coupled cluster module eT - Main program                                
!         Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017         
!                              
!                                             
!  :::::::::::::::::::::::::::::::::::::
!  -::- Modules used by the program -::-
!  :::::::::::::::::::::::::::::::::::::
!
!  IO module 
!
   use input_output
!
!  Method classes
!
   use hf_class
   use ccs_class
   use cc2_class
   use ccsd_class
!
   implicit none
!
!  Method class allocatable objects
!
   type(cc2), allocatable  :: cc2
   type(ccsd), allocatable :: ccsd
!
!  Wavefunction pointer
!
   class(hartree_fock), pointer :: wf
!
!  Unit identifier for input file eT.inp
!
   integer(i15) :: unit_identifier_input = -1
!
!  Open input file
!
   call generate_unit_identifier(unit_identifier_input)
   open(unit=unit_identifier_input, file='eT.inp', status='old', form='formatted')
   rewind(unit_identifier_input)
!
!  ::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Reading method section of input file -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::
!
!  Set method - Initialize object
!
!  call method_reader(STRING)
!
!  :::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Reading calculation section of input file -::- 
!  :::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Set calculation specifications
!
!  call calculation_reader(wf%tasks)
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Reading settings section of input file -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Set calculation settings settings
!
!  call settings_reader(wf%settings) 
!
!  Close input file
!
!  :::::::::::::::::::::::::::::
!  -::- Prepare output file -::- 
!  :::::::::::::::::::::::::::::
!
!  :::::::::::::::::::::::::
!  -::- Run Calculation -::- 
!  :::::::::::::::::::::::::
!
!  call wf%drv()
!
!  :::::::::::::::::::::::::::
!  -::- Close output file -::- 
!  :::::::::::::::::::::::::::
!
end program eT_program