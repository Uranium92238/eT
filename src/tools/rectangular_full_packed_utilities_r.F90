!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
module rectangular_full_packed_utilities_r
!
!!
!!    Rectangular full packed utilities (real)
!!
!!    Contains routines that perform various operations arrays
!!    in rectangular full packed (rfp) format.
!!    
!!    Rectangular full packed: 
!!    These routines only deal with the upper normal 
!!    version of the format. 
!!    The leading upper triangular matrix is transposed and tucked
!!    in with the trailing upper triangular matrix
!!
!!    x x x y y y    y y y      x x x y y y y    y y y y
!!      x x y y y    y y y        x x y y y y    y y y y
!!        x y y y -> y y y          x y y y y -> y y y y
!!          z z z    z z z            y y y y    y y y y
!!            z z    x z z              z z z    x z z z
!!              z    x x z                z z    x x z z
!!                   x x x                  z    x x x z
!!
!!    Note the difference between even and odd dimensions.
!!
!!    See: ACM Transactions on Mathematical Software, Vol. 37, No. 2, Article 18
!!         http://doi.acm.org/10.1145/1731022.1731028
!!
!
!!
!
   use parameters
!
contains
!
!
   pure subroutine rfp_index_to_full_index(x, y, dim, p, q)
!!
!!    RFP index to full index
!!    Written by Sarai D. Folkestad and Rolf H. Myhre, 2021
!!
      implicit none
!
      integer, intent(in) :: x, y, dim
      integer, intent(out) :: p, q
!
      if (x .le. y + dim/2) then
         p = x
         q = y + dim/2
      else
         p = y
         q = x - dim/2 - 1
      endif
!
   end subroutine rfp_index_to_full_index
!
!
   pure subroutine full_index_to_rfp_index(p, q, dim, x, y)
!!
!!    Full index to RFP index
!!    Written by Sarai D. Folkestad and Rolf H. Myhre, 2021
!!
      implicit none
!
      integer, intent(in) :: p, q, dim
      integer, intent(out) :: x, y
!
      integer :: pp, qq
!
      if (p .ge. q) then
!
         qq = p
         pp = q
!
      else
!
         qq = q
         pp = p
!
      endif
!
      if (qq .gt. dim/2) then
         y = qq - dim/2
         x = pp
      else
         y = pp
         x = qq + dim/2 + 1
      endif
!
   end subroutine full_index_to_rfp_index
!
!
   subroutine rfp_1234_to_2143(X_pqrs_rfp, Y_qpsr_rfp, dim_p, dim_q)
!!
!!    RFP 1234 to 2143
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Reorder array in rfp format from 1234 to 2143 order
!!
      implicit none
      integer, intent(in) :: dim_p, dim_q
!
      real(dp), dimension(dim_p*dim_q + mod(dim_p*dim_q+1,2), (dim_p*dim_q+1)/2), intent(in) :: X_pqrs_rfp
      real(dp), dimension(dim_p*dim_q + mod(dim_p*dim_q+1,2), (dim_p*dim_q+1)/2), intent(inout) ::  Y_qpsr_rfp
!
      integer :: dim, y, x, p, q, r, s, pq, rs, qp, sr, w, z
!
      dim = dim_p*dim_q
!
!$omp parallel do private(y, x, w, z, p, q, r, s, pq, qp, rs, sr)
      do y = 1, (dim+1)/2
!
         do x = 1, dim + mod(dim+1,2)
!
            call rfp_index_to_full_index(x, y, dim, pq, rs)
!
            p = mod(pq - 1,dim_p) + 1
            r = mod(rs - 1,dim_p) + 1
!
            q = (pq - 1)/dim_p + 1
            s = (rs - 1)/dim_p + 1
!
            qp = dim_q*(p-1) + q
            sr = dim_q*(r-1) + s
!
            call full_index_to_rfp_index(qp, sr, dim, w, z)

            Y_qpsr_rfp(w, z) = X_pqrs_rfp(x, y)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine rfp_1234_to_2143
!
!
   subroutine add_to_rfp(alpha, x, y, dim)
!!
!!    Add to RFP
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!       y = y + alpha x
!!
!!    where y is rfp and x is full symmetric
!!
      implicit none
!
      real(dp), intent(in) :: alpha
      integer, intent(in) :: dim
      real(dp), dimension(dim, dim), intent(in) :: x
      real(dp), dimension(dim + mod(dim+1,2), (dim+1)/2), intent(inout) :: y
!
      integer :: p, q, p_rfp, q_rfp
!
!$omp parallel do private (p, q, p_rfp, q_rfp)
      do p = 1, dim
         do q = 1, p
!
            call full_index_to_rfp_index(p, q, dim, p_rfp, q_rfp)
            y(p_rfp, q_rfp) = y(p_rfp, q_rfp) + alpha * x(p, q)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_to_rfp
!
!
end module rectangular_full_packed_utilities_r
