module excited_state_specs_class
!
!!
!!                       Excited State Specifictions class module                                 
!!            Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2017         
!!
!!    The class contains specifications for the excited state calculation, and necesary variables provided by user.
!! 
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
!
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the excited_state_specs class -::-
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   type :: excited_state_specs
!
!     Settings for the excited state calculation
!
      real(dp) :: energy_threshold   = 1.0D-06 ! Threshold for energies 
      real(dp) :: residual_threshold = 1.0D-06 ! Threshold for equation residuals
!
      integer(i15) :: max_iterations = 50 ! excited state maximum number of iterations 
!
      logical :: restart = .false.
!
      logical                                   :: user_specified_start_vector = .false. ! True if the user has provided startvectors
      integer(i15), dimension(:,:), allocatable :: start_vectors
!
!     Variables needed for excited state calculation
!
      integer(i15) :: n_singlet_states = 0
      integer(i15) :: n_triplet_states = 0
!
!     Logicals to determine wether we have Jacobi of Jacobi transpose transformation
!
      logical :: right = .true. ! Default
      logical :: left  = .false. 
!
!     Solution file name
!
      character(len=40) :: solution_file = 'right_valence' ! Should be 'right_valence', 'left_valence'
!                                                            'right_core' or 'left_core'.
!     Logical for printing excitation vectors 
!
      logical :: print_vectors = .true.
!
   end type excited_state_specs                                                                          
!
contains
!
end module excited_state_specs_class
