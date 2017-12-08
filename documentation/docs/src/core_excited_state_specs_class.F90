module core_excited_state_specs_class
!
!!
!!                  Core Excited State Specifictions class module                                 
!!            Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2017         
!!
!!    The class contains specifications for the core excited state calculation, and necesary variables provided by user.
!! 
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the core excited_state_specs class -::-
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   type :: core_excited_state_specs
!   
!     Variables necessary for core excitation calc
!
      integer(i15)                              :: n_equivalent_cores        ! Number of equivalent cores
      integer(i15), dimension(:,:), allocatable :: cores          ! Cores, given by order in MOLECULE.INP file
      integer(i15), dimension(:,:), allocatable :: index_core_mo  ! MO index/indices of core(s)
!
   end type core_excited_state_specs                                                                          
!
contains
!
end module core_excited_state_specs_class