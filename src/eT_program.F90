program eT_program
!
!
!                        eT coupled cluster program                                
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
!  Input file reader module 
!
   use input_reader
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
   class(hf), pointer :: wf => null()
!
!  Method string 
!
   character(len=25) :: method 
!
!  Unit identifier for input file eT.inp
!
   integer(i15) :: unit_input = -1
!
!  Open input file
!
   call generate_unit_identifier(unit_input)
   open(unit=unit_input, file='eT.inp', status='old', form='formatted')
   rewind(unit_input)
!
!  ::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Reading method section of input file -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::
!
!  Read the method from file, allocate the appropriate 
!  wavefunction object, and point to the main wavefunction 
!  wf to this object.
!
   call read_method(unit_input, method)
!
   if (method == 'CC2') then
!
      allocate(cc2)
      wf => cc2 
!
   elseif (method == 'CCSD') then
! 
      allocate(ccsd)
      wf => ccsd 
!
   else
!
      write(unit_output,*) 'Input error: method ', method, ' not recognized.'
      stop ! Terminate program
!
   endif
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