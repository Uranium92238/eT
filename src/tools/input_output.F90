module input_output
!
!!
!!    Input_output module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, April 2017
!!
!!    Handles input and output and contains:
!!
!!    generate_unit_identifier: Returns a free unit_identifier which can be assigned to a file that is opened.
!!    vec_print:                Prints vector in compound index order.
!!
!!    Note, June 2018, EFK. THIS MODULE SHOULD BE REMOVED ONCE ALL CALLS TO
!!    GENERATE UNIT IDENTIFIER ARE GONE. REMEMBER TO MOVE "UNIT_OUTPUT"
!!    TO SOME OTHER POSITION WHEN DELETING THE MODULE.
!!
!
   use types
!
   implicit none
!
   integer(i15) :: unit_output = 0
   integer, private :: n_files = 0
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
   subroutine vec_print(vec,dim_1,dim_2)
!!
!!    Vector print
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, March 2017
!!
!!    A vector is printed with a compound index (p q) of dimension (dim_1 x dim_2)
!!
      implicit none
!
      integer(i15) :: p = 0, q = 0, pq = 0
!
      integer(i15), intent(in) :: dim_1,dim_2
      real(dp), dimension(dim_1, dim_2), intent(in) :: vec
!
      do q = 1, dim_2
         do p = 1, dim_1
!
            write(unit_output,*) p, q, vec(p,q)
!
         enddo
      enddo
!
   end subroutine vec_print
!
!
   subroutine vec_print_nonzero_elm(vec,dim_1,dim_2)
!!
!!    Vector print
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, March 2017
!!
!!    A vector is printed with a compound index (p q) of dimension (dim_1 x dim_2)
!!
      implicit none
!
      integer(i15) :: p = 0, q = 0, pq = 0
!
      integer(i15), intent(in) :: dim_1,dim_2
      real(dp), dimension(dim_1, dim_2), intent(in) :: vec
!
      do q = 1, dim_2
         do p = 1, dim_1
!
            if (vec(p, q) .gt. 1.0D-03) then
               write(unit_output,*) p, q, vec(p,q)
            endif
!
         enddo
      enddo
!
   end subroutine vec_print_nonzero_elm
end module input_output
