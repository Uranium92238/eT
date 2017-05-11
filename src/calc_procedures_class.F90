module calc_procedures_class
!
!
!                         Calculation procedures class module                                 
!              Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017         
!
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
!
!  :::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the calc_procedures class -::-
!  :::::::::::::::::::::::::::::::::::::::::::::::
!
   type :: calc_procedures 
!
      logical :: do_ground_state   = .false.
      logical :: do_excited_state  = .false. 
      logical :: do_properties     = .false.
!
      integer(i15) :: n_singlet_states = 0
      integer(i15) :: n_triplet_states = 0
!
   end type calc_procedures                                                                            
!
!
contains
!
end module calc_procedures_class