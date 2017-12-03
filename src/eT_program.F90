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
   use mp2_class
   use ccs_class
   use cc2_class
   use ccsd_class
   use cc3_class
   use ccsdpt_class
   use mlcc2_class
   use mlccsd_class
!
   implicit none
!
!  Method class allocatable objects
!
   type(mp2),    allocatable, target :: mp2_wf
   type(ccs),    allocatable, target :: ccs_wf 
   type(cc2),    allocatable, target :: cc2_wf
   type(ccsd),   allocatable, target :: ccsd_wf
   type(cc3),    allocatable, target :: cc3_wf
   type(ccsdpt), allocatable, target :: ccsdpt_wf
!
!  Multi-level methods
!
   type(mlcc2), allocatable, target  :: mlcc2_wf
   type(mlccsd), allocatable, target :: mlccsd_wf
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
   character(len=13) :: orbital_level
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
   if (trim(method) == 'MP2') then 
!
      allocate(mp2_wf)
      wf => mp2_wf
!
   elseif (trim(method) == 'CCS') then 
!
      allocate(ccs_wf)
      wf => ccs_wf
!
   elseif (trim(method) == 'CCSD') then
!
      allocate(ccsd_wf)
      wf => ccsd_wf
!
   elseif (trim(method) == 'CC2') then
!
      allocate(cc2_wf)
      wf => cc2_wf
!
   elseif (trim(method) == 'CC3') then 
!
      allocate (cc3_wf)
      wf => cc3_wf
!
   elseif (trim(method) == 'CCSD(T)') then
!
      allocate(ccsdpt_wf)
      wf => ccsdpt_wf
!
   elseif (trim(method) == 'MLCC2') then
!
      allocate(mlcc2_wf)
      wf => mlcc2_wf
!
   elseif (trim(method) == 'MLCCSD') then
!
      allocate(mlccsd_wf)
      wf => mlccsd_wf
!
   else
      write(unit_output,*) 'Method ', trim(method), ' not recognized.'
      flush(unit_output)
      stop
!
   endif
!
! Close input file
!
   close(unit_input)
!
! :::::::::::::::::::::::::
! -::- Run Calculation -::- 
! :::::::::::::::::::::::::
!
   call wf%init
   call wf%drv
!
! :::::::::::::::::::::::::::
! -::- Close output file -::- 
! :::::::::::::::::::::::::::
!
  close(unit_output)
!
end program eT_program
