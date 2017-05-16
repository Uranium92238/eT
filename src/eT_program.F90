program eT_program
!
!!
!!                        eT - a coupled cluster program                                
!!         Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017         
!!                              
!                                             
!  :::::::::::::::::::::::::::::::::::::
!  -::- Modules used by the program -::-
!  :::::::::::::::::::::::::::::::::::::
!
   use input_output ! IO module 
   use input_reader ! Input reader routines
!
!  Wavefunction classes
!
   use hf_class
   use ccs_class
   use ccsd_class
!
   implicit none
!
!  Method class allocatable objects
!
   type(ccs),  allocatable, target :: ccs_wf 
   type(ccsd), allocatable, target :: ccsd_wf
!
!  Wavefunction pointer
!
   class(hf), pointer :: wf => null()
!
!  Method string 
!
   character(len=40) :: method 
!
!  Unit identifier for input file eT.inp
!
   integer(i15) :: unit_input = -1
!
!  ::::::::::::::::::::::::::::::::::::::::::::::
!  -::-     Set up memory controller         -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::
!
   call work_init
!
!  ::::::::::::::::::::::::::::::::::::::::::::::
!  -::-  Create & open the main output file  -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::
!
   call generate_unit_identifier(unit_output)
   open(unit=unit_output,file='eT.out',status='unknown',form='formatted')
   rewind(unit_output)
!
!  ::::::::::::::::::::::::::::::::::::::::::::::
!  -::-         Print program banner         -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::
!
   write(unit_output,'(//t18,a)')  'eT - a coupled cluster program'
   write(unit_output,'(t15,a//)') 'S. D. Folkestad, E. F. Kjønstad, 2017'
!
!  ::::::::::::::::::::::::::::::::::::::::::::::
!  -::-    Print banner for input section    -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::
!
   write(unit_output,'(t3,a/)')   ':: Input reader'
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
   call method_reader(unit_input, method)
!
   write(unit_output,'(t3,a,a,a)') 'Our wavefunction is of type ',trim(method),'.'
   flush(unit_output)
!
   if (trim(method) == 'CCS') then 
!
      allocate(ccs_wf)
      wf => ccs_wf
!
   elseif (trim(method) == 'CCSD') then
!
      allocate(ccsd_wf)
      wf => ccsd_wf
!
   else
!
      write(unit_output,*) 'Method ', trim(method), ' not recognized.'
      flush(unit_output)
      stop
!
   endif
!
!  :::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Reading calculation section of input file -::- 
!  :::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Set calculation specifications
!
   call calculation_reader(unit_input, wf%tasks)
!
   if (wf%tasks%ground_state)  write(unit_output,'(t3,a)')  'Ground state calculation requested.'
   if (wf%tasks%excited_state) write(unit_output,'(t3,a)')  'Excited state calculation requested.' ! Dummy as of now 
   if (wf%tasks%properties)    write(unit_output,'(t3,a)')  'Properties calculation requested.'    ! Dummy as of now
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Reading settings section of input file -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Set the calculation settings of the wavefunction
!
   call settings_reader(unit_input, wf%settings) 
!
   write(unit_output,'(/t3,a/)')         'Settings for this calculation:'
!
   write(unit_output,'(t6,a25,e14.2)')    'Energy threshold:',         wf%settings%energy_threshold
   write(unit_output,'(t6,a25,e14.2)')    'Amplitude eqs. threshold:', wf%settings%ampeqs_threshold
   write(unit_output,'(t6,a25,i14/)')     'Memory:',                   mem
!
!  Close input file
!
   close(unit_input)
!
!  :::::::::::::::::::::::::
!  -::- Run Calculation -::- 
!  :::::::::::::::::::::::::
!
   call wf%init
   call wf%drv
!
!  :::::::::::::::::::::::::::
!  -::- Close output file -::- 
!  :::::::::::::::::::::::::::
!
   close(unit_output)
!
end program eT_program
