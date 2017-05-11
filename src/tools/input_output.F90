module input_output
!
!!
!!    Input_output module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, April 2017
!!
!!    Module handles input and output and contains:
!! 
!!    generate_unit_identifier: Returns a free unit_identifier which can be assigned to a file that is opened. 
!!    vec_print:                Prints vector in compound index order.  
!!            
!
   use types
!
   implicit none
!
   integer(i15) :: unit_output = 0 
   integer, private :: n_files = 200 ! To not overwrite the DALTON.OUT identifier (for debug against old code)
!
contains
!
!
   subroutine generate_unit_identifier(unit_identifier)
!!
!!    Generate unit identifier
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, April 2017
!!
!!    Returns a valid unit identifier for opening file. 
!!
      implicit none
!
      integer(i15) :: unit_identifier
!
      n_files = n_files + 1
      unit_identifier = n_files
!
   end subroutine generate_unit_identifier
!
!
end module input_output
