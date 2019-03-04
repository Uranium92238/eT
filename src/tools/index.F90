!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
module index
!!
!!    Index module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2018
!!
!!    This module contains routines to calculate compound indices, such as pq, or pqr,
!!    along with routines to calculate packed indices for symmetric matrices, and also,
!!    inversion of indices (from compound indices to individual indices).
!!
!
   use kinds
   use parameters
   use file_class
!
   implicit none
!
contains
!
!
   integer function index_packed(i,j)
!!
!!    Packed index
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Returns a packed index (used for symmetric matrices)
!!
      implicit none
!
      integer, intent(in) :: i,j
!
      index_packed = (max(i,j)*(max(i,j)-3)/2) + i + j
!
   end function index_packed
!
!
   integer function index_three(p,q,r,dim_p,dim_q)
!!
!!    Three index compound
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Returns the compound index (pqr)
!!
      implicit none
!
      integer, intent(in) :: p, q, r, dim_p, dim_q
!
      index_three = dim_p*(dim_q*(r-1)+q-1)+p
!
   end function index_three
!
!
   integer function index_two(p, q, dim_p)
!!
!!    Two index compound
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Returns the compound index (pq)
!!
      implicit none
!
      integer, intent(in) :: p, q, dim_p
!
      index_two = dim_p*(q-1)+p
!
   end function index_two
!
!
   integer function packed_size(N)
!!
!!    Packed size
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Returns size of packed symmetric matrices
!!    of dimension N x N (triangular elements)
!!
      implicit none
!
      integer, intent(in) :: N
!
      packed_size = N*(N+1)/2
!
   end function packed_size
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
      integer, intent(in) :: pq
!
      integer :: p, q
      integer, intent(in) :: dim_p, dim_q
!
      integer :: I = 0
!
!     Since dim_p*q >= pq > dim_p*(q-1) by construction, we can determine q
!
      do I = 1, dim_q
!
         if  (pq .gt. dim_p*(I-1) .and. pq .le. dim_p*I) then
!
            q = I
!
         endif
!
      enddo
!
!     From q, it is straight-forward to determine p
!
      p = pq - dim_p*(q-1)
!
!     Sanity test
!
      if (index_two(p, q, dim_p) .ne. pq) then
!
         write(output%unit,'(/t3,a)') 'Error: inversion of compound index unsuccessful'
         stop
!
      endif
!
   end subroutine invert_compound_index
!
!
   subroutine invert_packed_index(pq, p, q, dim)
!!
!!    Invert packed index
!!    Written by Eirik F. Kjønstad, June 2018
!!
!!    Returns the indices (p, q) in the upper triangular part of the symmetric matrix.
!!
      implicit none
!
      integer :: pq
!
      integer :: p, q
!
      integer :: dim
!
      integer :: I = 0, J = 0
!
!     Loop through upper triangular part until the index matches
!
      do I = 1, dim
         do J = I, dim
!
            if (index_packed(I, J) .eq. pq) then
!
               p = I
               q = J
!
            endif
!
         enddo
      enddo
!
!     Sanity check
!
      if (index_packed(p, q) .ne. pq) then
!
         write(output%unit,'(/t6,a)') 'Error: inversion of packed index unsuccessful'
         stop
!
      endif
!
   end subroutine invert_packed_index
!
!
end module index
