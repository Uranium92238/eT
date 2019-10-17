module index_invert
!
!!
!!    Index invert module
!!
!!    Routines to invert compound and packed compound indices.
!!
!
   use kinds
   use parameters
!
   implicit none
!
contains
!
!
   subroutine invert_compound_index(pq, p, q, dim_p, dim_q)
!!
!!    Invert compound index
!!    Written by Eirik F. Kjønstad, May 2018
!!
!!    Given the compound index pq = dim_p*(q-1) + p, this routine
!!    determines p and q.
!!
      implicit none
!
      integer, intent(in) :: pq, dim_p, dim_q
      integer, intent(out) :: p, q
      integer :: I
!
!     Since dim_p*q >= pq > dim_p*(q-1) by construction, we can determine q
!
      do I = 1, dim_q
         if  (pq .gt. dim_p*(I-1) .and. pq .le. dim_p*I) then
            q = I
            exit
         endif
      enddo
!
!     From q, it is straight-forward to determine p
!
      p = pq - dim_p*(q-1)
!
   end subroutine invert_compound_index
!
!
   subroutine invert_packed_index(pq, p, q, dim_)
!!
!!    Invert packed index
!!    Written by Eirik F. Kjønstad, June 2018
!!
!!    Returns the indices (p, q) in the upper triangular part of the symmetric matrix.
!!
      implicit none
!
      integer, intent(in)  :: pq, dim_
      integer, intent(out) :: p, q
      integer :: I, J
      logical :: found
!
!     Loop through upper triangular part until the index matches
!
      found = .false.
      do I = 1, dim_
         do J = I, dim_
            if (((J*J-3*J)/2 + J + I) .eq. pq) then
               p = I
               q = J
               found = .true.
               exit
            endif
         enddo
         if (found) exit
      enddo
!
   end subroutine invert_packed_index
!
!
end module index_invert
