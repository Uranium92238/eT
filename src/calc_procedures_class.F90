module calc_procedures_class
!
!!
!!                       Calculation procedures class module                                 
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
   type :: calc_procedures 
!
      logical :: ground_state       = .false.
      logical :: excited_state      = .false. 
      logical :: core_excited_state = .false. 
      logical :: properties         = .false.
!
!     For excited states (core and valence)
!
      integer(i15) :: n_singlet_states = 0
      integer(i15) :: n_triplet_states = 0
!
!     For core excited states 
!
      integer(i15) :: n_cores ! number of equivalent cores
      integer(i15), dimension(:,:), allocatable :: cores
      integer(i15), dimension(:,:), allocatable :: index_core_mo
!
   end type calc_procedures                                                                            
!
contains
!
end module calc_procedures_class