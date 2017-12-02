module calc_tasks_class
!
!!
!!                       Calculation tasks class module                                 
!!            Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017         
!!
!!    The procedures class is used for two purposes. Wavefunctions possess 
!!    "implemented" and "tasks" objects, which are both instances of the class 
!!    calculation procedures. If one of the tasks requested is not implemented
!!    for the wavefunction in question, the program stops, letting the user know
!!    that the requested functionality is not a part of eT (yet).
!!
!!    Tasks are determined during the read of the input (see the main program).
!!    Implemented are set by the wavefunction's initialization routine.
!! 
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
!
!  :::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the calc_procedures class -::-
!  :::::::::::::::::::::::::::::::::::::::::::::::::
!
   type :: calc_tasks
!
      logical :: ground_state       = .false.
      logical :: excited_state      = .false. 
      logical :: core_excited_state = .false. 
      logical :: ionized_state      = .false. 
      logical :: core_ionized_state = .false. 
      logical :: properties         = .false.
!
      character(len=40) :: current = 'ground_state'
!
   end type calc_tasks                                                                           
!
contains
!
end module calc_tasks_class