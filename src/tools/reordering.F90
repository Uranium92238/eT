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
module reordering
!
!!
!!    Reordering module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2018
!!
!!    Here are all routines that are used to re-sort arrays and (e.g., sort_123_to_132)
!!    performing re-sort-addition (add_132_to_123), along with routines that both squareup
!!    (e.g., t2am -> t_ai_bj) and re-sort (t_ai_bj -> t_aj_bi) in a single step; also included
!!    are routines that squareup (t_aibj -> t_ai_bj) and packin (t_ai_bj -> t_aibj) and
!!    related routines.
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
!     -::- Two-index re-sort and re-sort-add routines -::-
!     ----------------------------------------------------
!
   subroutine sort_12_to_21(x_p_q, x_q_p, dim_p, dim_q)
!!
!!    Sort 12 to 21
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_p_q to x_q_p (i.e., 12 to 21).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_sp_rq is assumed allocated as dim_s*dim_p x dim_r*dim_q.
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q
!
      real(dp), dimension(dim_p, dim_q), intent(in) :: x_p_q
      real(dp), dimension(dim_q, dim_p) :: x_q_p
!
      integer :: p, q
!
!$omp parallel do schedule(static) private(p,q)
      do p = 1, dim_p
         do q = 1, dim_q
!
            x_q_p(q, p) = x_p_q(p, q)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_12_to_21
!
!
   subroutine add_21_to_12(scalar, x, y_p_q, dim_p, dim_q)
!!
!!    Add 21 to 12
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_p_q(p,q) = y_p_q(p,q) + scalar*x(q,p)
!!
!!    The unordered array y_p_q is assumed allocated as dim_p x dim_q,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q
!
      real(dp), dimension(dim_p, dim_q) :: y_p_q
      real(dp), dimension(dim_q, dim_p), intent(in) :: x
!
      integer :: p, q
!
      do q = 1, dim_q
         do p = 1, dim_p
!
               y_p_q(p,q) = y_p_q(p,q) + scalar*x(q,p)
!
         enddo
      enddo
!
   end subroutine add_21_to_12
!
!
   subroutine symmetric_sum(x, dim_)
!!
!!    Symmetric sum
!!    Written by Eirik F. Kjønstad, Dec 2017
!!
!!    Performs the action
!!
!!       x(p,q) = x(p,q) + x(q,p)
!!
!!    without making a separate copy of x.
!!
!!    Note: the effect is equal to the routine "add_21_to_12" with scalar = 1,
!!          where a second copy of x is not needed.
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_, dim_), intent(inout) :: x
!
      integer :: p, q
!
!     Overwrite the lower triangular part of the matrix
!
!$omp parallel do private(p, q)
      do q = 1, dim_
         do p = q, dim_
!
            x(p,q) = x(p,q) + x(q,p)
!
         enddo
      enddo
!$omp end parallel do
!
!     Copy the lower triangular part to the upper triangular part
!
!$omp parallel do private(p, q)
      do p = 1, dim_
         do q = p + 1, dim_
!
            x(p,q) = x(q,p)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine symmetric_sum
!
!
!     -::- Three-index re-sort and re-sort-add routines -::-
!     ------------------------------------------------------
!
   subroutine sort_123_to_312(x_pqr, x_rpq, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 312
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       x_rpq(r,p,q) = x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r), intent(in)  :: x_pqr
      real(dp), dimension(dim_r, dim_p, dim_q), intent(out) :: x_rpq
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do q = 1, dim_q
         do p = 1, dim_p
            do r = 1, dim_r
!
               x_rpq(r,p,q) = x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_312
!
!
   subroutine sort_123_to_231(x_pqr, x_qrp, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 312
!!    Written by Rolf H. Myhre, April 2019
!!
!!    Performs:
!!
!!       x_qrp(q,r,p) = x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r), intent(in)  :: x_pqr
      real(dp), dimension(dim_q, dim_r, dim_p), intent(out) :: x_qrp
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do p = 1, dim_p
         do r = 1, dim_r
            do q = 1, dim_q
!
               x_qrp(q,r,p) = x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_231
!
!
   subroutine sort_123_to_312_and_add(x_pqr, x_rpq, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 312 and add
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, January 2019
!!
!!    Performs:
!!
!!       x_rpq(r,p,q) = x_rpq(r,p,q) + x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r), intent(in)    :: x_pqr
      real(dp), dimension(dim_r, dim_p, dim_q), intent(inout) :: x_rpq
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do q = 1, dim_q
         do p = 1, dim_p
            do r = 1, dim_r
!
               x_rpq(r,p,q) = x_rpq(r,p,q) + x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_312_and_add
!
!
   subroutine sort_123_to_231_and_add(x_pqr, x_qrp, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 231 and add
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, January 2019
!!
!!    Performs:
!!
!!       x_qrp(q,r,p) = x_qrp(q,r,p) + x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r), intent(in)    :: x_pqr
      real(dp), dimension(dim_q, dim_r, dim_p), intent(inout) :: x_qrp
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do p = 1, dim_p
         do r = 1, dim_r
            do q = 1, dim_q
!
               x_qrp(q,r,p) = x_qrp(q,r,p) + x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_231_and_add
!
!
   subroutine sort_123_to_321(x_pqr, x_rqp, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 321
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       x_rqp(r,q,p) = x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r), intent(in)    :: x_pqr
      real(dp), dimension(dim_r, dim_q, dim_p), intent(inout) :: x_rqp
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do p = 1, dim_p
         do q = 1, dim_q
            do r = 1, dim_r
!
               x_rqp(r,q,p) = x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_321
!
!
   subroutine construct_123_minus_321(x_pqr, y_pqr, dim_)
!!
!!    Construct 123 minus 321
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, January 2019
!!
!!    Performs:
!!
!!       y_pqr(p,q,r) = x_pqr(p,q,r) - x_pqr(r,q,p)
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_,dim_,dim_), intent(in)  :: x_pqr
      real(dp), dimension(dim_,dim_,dim_), intent(out) :: y_pqr
!
      integer :: p,q,r
!
!$omp parallel do schedule(static) private(p,q,r)
      do r = 1, dim_
         do q = 1, dim_
           do p = 1, dim_
!
               y_pqr(p,q,r) = x_pqr(p,q,r) - x_pqr(r,q,p)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_123_minus_321
!
!
   subroutine construct_123_minus_132(x_pqr, y_pqr, dim_)
!!
!!    Construct 123 minus 132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, January 2019
!!
!!    Performs:
!!
!!       y_pqr(p,q,r) = x_pqr(p,q,r) - x_pqr(p,r,q)
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_,dim_,dim_), intent(in)  :: x_pqr
      real(dp), dimension(dim_,dim_,dim_), intent(out) :: y_pqr
!
      integer :: p,q,r
!
!$omp parallel do schedule(static) private(p,q,r)
      do r = 1, dim_
         do q = 1, dim_
           do p = 1, dim_
!
               y_pqr(p,q,r) = x_pqr(p,q,r) - x_pqr(p,r,q)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_123_minus_132
!
!
   subroutine construct_123_min_132_min_321(x_pqr, y_pqr, dim_)
!!
!!    Construct 123 minus 132 minus 321
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, January 2019
!!
!!    Performs:
!!
!!       y_pqr(p,q,r) = 2*x_pqr(p,q,r) - x_pqr(p,r,q) - x_pqr(r,q,p)
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_,dim_,dim_), intent(in)  :: x_pqr
      real(dp), dimension(dim_,dim_,dim_), intent(out) :: y_pqr
!
      integer :: p,q,r
!
!$omp parallel do schedule(static) private(p,q,r)
      do r = 1, dim_
         do q = 1, dim_
            do p = 1, dim_
!
               y_pqr(p,q,r) = 2*x_pqr(p,q,r) - x_pqr(p,r,q) - x_pqr(r,q,p)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_123_min_132_min_321
!
!
   subroutine construct_132_minus_312(x_pqr, y_pqr, dim_)
!!
!!    Construct 132 minus 312
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, January 2019
!!
!!    Performs:
!!
!!       y_pqr(p,q,r) = x_pqr(p,r,q) - x_pqr(r,p,q)
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_,dim_,dim_), intent(in)  :: x_pqr
      real(dp), dimension(dim_,dim_,dim_), intent(out) :: y_pqr
!
      integer :: p,q,r
!
!$omp parallel do schedule(static) private(p,q,r)
      do r = 1, dim_
         do q = 1, dim_
            do p = 1, dim_
!
               y_pqr(p,q,r) = x_pqr(p,r,q) - x_pqr(r,p,q)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_132_minus_312
!
!
   subroutine construct_132_min_123_min_312(x_pqr, y_pqr, dim_)
!!
!!    Construct 132 minus 123 minus 312
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, January 2019
!!
!!    Performs:
!!
!!       y_pqr(p,q,r) = 2*x_pqr(p,r,q)  - x_pqr(p,q,r) - x_pqr(r,p,q)
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_,dim_,dim_), intent(in)  :: x_pqr
      real(dp), dimension(dim_,dim_,dim_), intent(out) :: y_pqr
!
      integer :: p,q,r
!
!$omp parallel do schedule(static) private(p,q,r)
      do r = 1, dim_
         do q = 1, dim_
            do p = 1, dim_
!
               y_pqr(p,q,r) = 2*x_pqr(p,r,q) - x_pqr(p,q,r) - x_pqr(r,p,q)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_132_min_123_min_312
!
!
   subroutine construct_321_minus_231(x_pqr, y_pqr, dim_)
!!
!!    Construct 321 minus 231
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, January 2019
!!
!!    Performs:
!!
!!       y_pqr(p,q,r) = x_pqr(r,q,p) - x_pqr(q,r,p)
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_,dim_,dim_), intent(in)  :: x_pqr
      real(dp), dimension(dim_,dim_,dim_), intent(out) :: y_pqr
!
      integer :: p,q,r
!
!$omp parallel do schedule(static) private(p,q,r)
      do r = 1, dim_
         do q = 1, dim_
            do p = 1, dim_
!
               y_pqr(p,q,r) = x_pqr(r,q,p) - x_pqr(q,r,p)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_321_minus_231
!
!
   subroutine construct_321_min_231_min_123(x_pqr, y_pqr, dim_)
!!
!!    Construct 321 minus 231 minus 123
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, January 2019
!!
!!    Performs:
!!
!!       y_pqr(p,q,r) = 2*x_pqr(r,q,p)  - x_pqr(q,r,p) - x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_,dim_,dim_), intent(in)  :: x_pqr
      real(dp), dimension(dim_,dim_,dim_), intent(out) :: y_pqr
!
      integer :: p,q,r
!
!$omp parallel do schedule(static) private(p,q,r)
      do r = 1, dim_
         do q = 1, dim_
            do p = 1, dim_
!
               y_pqr(p,q,r) = 2*x_pqr(r,q,p) - x_pqr(q,r,p) - x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_321_min_231_min_123
!
!
   subroutine construct_213_minus_231(x_pqr, y_pqr, dim_)
!!
!!    Construct 213 minus 231
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, January 2019
!!
!!    Performs:
!!
!!       y_pqr(p,q,r) = x_pqr(p,q,r) - x_pqr(r,q,p)
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_,dim_,dim_), intent(in)  :: x_pqr
      real(dp), dimension(dim_,dim_,dim_), intent(out) :: y_pqr
!
      integer :: p,q,r
!
!$omp parallel do schedule(static) private(p,q,r)
      do r = 1, dim_
         do q = 1, dim_
            do p = 1, dim_
!
               y_pqr(p,q,r) = x_pqr(q,p,r) - x_pqr(q,r,p)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_213_minus_231
!
!
   subroutine construct_213_minus_312(x_pqr, y_pqr, dim_)
!!
!!    Construct 213 minus 312
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, January 2019
!!
!!    Performs:
!!
!!       y_pqr(p,q,r) = x_pqr(q,p,r) - x_pqr(r,p,q)
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_,dim_,dim_), intent(in)  :: x_pqr
      real(dp), dimension(dim_,dim_,dim_), intent(out) :: y_pqr
!
      integer :: p,q,r
!
!$omp parallel do schedule(static) private(p,q,r)
      do r = 1, dim_
         do q = 1, dim_
            do p = 1, dim_
!
               y_pqr(p,q,r) = x_pqr(q,p,r) - x_pqr(r,p,q)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_213_minus_312
!
!
   subroutine sort_123_to_321_and_add(x_pqr, x_rqp, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 321 and add
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, January 2019
!!
!!    Performs:
!!
!!       x_rqp(r,q,p) = x_rqp(r,q,p) + x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r), intent(in)    :: x_pqr
      real(dp), dimension(dim_r, dim_q, dim_p), intent(inout) :: x_rqp
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r, q, p)
      do p = 1, dim_p
         do q = 1, dim_q
            do r = 1, dim_r
!
              x_rqp(r,q,p) = x_rqp(r,q,p) + x_pqr(p,q,r)
!
           enddo
        enddo
     enddo
!$omp end parallel do
!
   end subroutine sort_123_to_321_and_add
!
!
   subroutine sort_123_to_132(x_pqr, x_prq, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       x_prq(p,r,q) = x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r), intent(in)  :: x_pqr
      real(dp), dimension(dim_p, dim_r, dim_q), intent(out) :: x_prq
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do q = 1, dim_q
         do r = 1, dim_r
            do p = 1, dim_p
!
               x_prq(p,r,q) = x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_132
!
!
   subroutine sort_123_to_132_and_add(x_pqr, x_prq, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 132 and add
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, January 2019
!!
!!    Performs:
!!
!!       x_prq(p,r,q) = x_prq(p,r,q) + x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r), intent(in)    :: x_pqr
      real(dp), dimension(dim_p, dim_r, dim_q), intent(inout) :: x_prq
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do q = 1, dim_q
         do r = 1, dim_r
            do p = 1, dim_p
!
               x_prq(p,r,q) = x_prq(p,r,q) + x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_132_and_add
!
!
   subroutine sort_123_to_213(x_pqr, x_qpr, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       x_qpr(q,p,r) = x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r), intent(in)  :: x_pqr
      real(dp), dimension(dim_q, dim_p, dim_r), intent(out) :: x_qpr
!
      integer :: p, q, r
!
!$omp parallel do schedule(static) private(p, q, r)
      do r = 1, dim_r
         do p = 1, dim_p
            do q = 1, dim_q
!
               x_qpr(q,p,r) = x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_213
!
!
   subroutine sort_123_to_213_and_add(x_pqr, x_qpr, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 132 and add
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, January 2019
!!
!!    Performs:
!!
!!       x_qpr(q,p,r) = x_qpr(q,p,r) + x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p,dim_q,dim_r), intent(in)    :: x_pqr
      real(dp), dimension(dim_q,dim_p,dim_r), intent(inout) :: x_qpr
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do r = 1, dim_r
         do p = 1, dim_p
            do q = 1, dim_q
!
               x_qpr(q,p,r) = x_qpr(q,p,r) +  x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_213_and_add
!
!
   subroutine add_213_to_123(scalar, x_qpr, y_pqr, dim_p, dim_q, dim_r)
!!
!!    Add 213 to 123
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqr(pqr) = y_pqr(pqr) + scalar*x(qpr)
!!
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r) :: y_pqr
      real(dp), dimension(dim_q, dim_p, dim_r), intent(in) :: x_qpr
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do r = 1, dim_r
         do q = 1, dim_q
            do p = 1, dim_p
!
               y_pqr(p,q,r) = y_pqr(p,q,r) + scalar*x_qpr(q,p,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_213_to_123
!
!
   subroutine add_132_to_123(scalar, x_prq, y_pqr, dim_p, dim_q, dim_r)
!!
!!    Add 132 to 123
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqr(p,q,r) = y_pqr(p,q,r) + scalar * x_prq(p,r,q)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r) :: y_pqr
      real(dp), dimension(dim_p, dim_r, dim_q), intent(in) :: x_prq
!
      integer :: p, q, r
!
!$omp parallel do schedule(static) private(r,q,p)
      do r = 1, dim_r
         do q = 1, dim_q
            do p = 1, dim_p
!
               y_pqr(p,q,r) = y_pqr(p,q,r) + scalar*x_prq(p,r,q)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_132_to_123
!
!
!     -::- Four-index re-sort and re-sort-add routines -::-
!     -----------------------------------------------------
!
   subroutine add_3124_to_1234(scalar, x_rpqs, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 3124 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_rpqs(r,p,q,s)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_r, dim_p, dim_q, dim_s), intent(in) :: x_rpqs
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_rpqs(r,p,q,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_3124_to_1234
!
!
   subroutine sort_1234_to_3412(x_pqrs, x_rspq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3412
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_rspq (i.e., 1234 to 3412).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_r, dim_s, dim_p, dim_q), intent(out) :: x_rspq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do p = 1, dim_p
            do s = 1, dim_s
               do r = 1, dim_r
!
                  x_rspq(r, s, p, q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_3412
!
!
   subroutine sort_1234_to_4132(x_pqrs, x_sprq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_sprq (i.e., 1234 to 4132).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_s, dim_p, dim_r, dim_q), intent(out) :: x_sprq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do r = 1, dim_r
            do p = 1, dim_p
               do s = 1, dim_s
!
                  x_sprq(s,p,r,q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_4132
!
   subroutine squareup_and_sort_1234_to_4132(x_pqrs, x_sprq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Square up and sort 1234 to 4132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre
!!    and Andreas Skeidsvoll, 2018
!!
!!    Reorders the array x_pq_rs to x_sp_rq (i.e., 1234 to 4132).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_sp_rq is assumed allocated as dim_s*dim_p x dim_r*dim_q.
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*(dim_r*dim_s))/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_s, dim_p, dim_r, dim_q), intent(out)       :: x_sprq
!
      integer :: p, q, r, s, rs, pq, pqrs
!
!$omp parallel do schedule(static) private(s,r,rs,q,p,pq,pqrs)
      do q = 1, dim_q
         do r = 1, dim_r
            do p = 1, dim_p
!
               pq = dim_p*(q-1) + p
!
               do s = 1, dim_s
!
                  rs = dim_r*(s-1) + r
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_sprq(s,p,r,q) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_4132
!
!
   subroutine sort_1234_to_4123(x_pqrs, x_spqr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4132
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and Andreas Skeidsvoll, 2018
!!
!!    Reorders the array x_pqrs to x_spqr (i.e., 1234 to 4123).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_s, dim_p, dim_q, dim_r), intent(out) :: x_spqr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do r = 1, dim_r
            do p = 1, dim_p
               do s = 1, dim_s
!
                  x_spqr(s,p,q,r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_4123
!
!
   subroutine squareup_and_sort_1234_to_4123(x_pqrs, x_spqr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Square up and sort 1234 to 4132
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and Andreas Skeidsvoll, 2018
!!
!!    Reorders and unpacks the array x_pqrs to x_spqr (i.e., 1234 to 4123).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*(dim_r*dim_s))/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_s, dim_p, dim_q, dim_r), intent(out)       :: x_spqr
!
      integer :: p, q, r, s, rs, pq, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do r = 1, dim_r
         do q = 1, dim_q
            do p = 1, dim_p
!
               pq = dim_p*(q-1) + p
!
               do s = 1, dim_s
!
                  rs = dim_r*(s-1) + r
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_spqr(s,p,q,r) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_4123
!
!
   subroutine sort_1234_to_3124(x_pqrs, x_rpqs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3124
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_rpqs (i.e., 1234 to 3124).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_r, dim_p, dim_q, dim_s), intent(out) :: x_rpqs
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do q = 1, dim_q
            do p = 1, dim_p
               do r = 1, dim_r
!
                  x_rpqs(r,p,q,s) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_3124
!
!
   subroutine sort_1234_to_3142(x_pqrs, x_rpsq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3142
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and
!!    Sarai D. Folkestad, 2018
!!
!!    Reorders the array x_pqrs to x_rpsq (i.e., 1234 to 3142).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_r, dim_p, dim_s, dim_q), intent(out) :: x_rpsq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do s = 1, dim_s
            do p = 1, dim_p
               do r = 1, dim_r
!
                  x_rpsq(r,p,s,q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_3142
!
!
   subroutine squareup_and_sort_1234_to_2413(x_pqrs, x_qspr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 2413
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Square up and reorder the array x_pqrs to x_qspr (i.e., 1234 to 2413).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*dim_r*dim_s)/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_q, dim_s, dim_p, dim_r), intent(out)     :: x_qspr
!
      integer :: p, q, r, s, rs, pq, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do r = 1, dim_r
         do p = 1, dim_p
            do s = 1, dim_s
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
                  pq = dim_p*(q-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_qspr(q,s,p,r) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_2413
!
!
   subroutine squareup_and_sort_1234_to_2341(x_pqrs, x_qrsp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 2341
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!
!!    Unpack and reorder the array x_pqrs to x_qrsp (i.e., 1234 to 2413).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*(dim_r*dim_s))/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_q, dim_r, dim_s, dim_p), intent(out)       :: x_qrsp
!
      integer :: p, q, r, s, rs, pq, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do p = 1, dim_p
         do s = 1, dim_s
            do r = 1, dim_r
!
               rs = dim_r*(s-1) + r
!
               do q = 1, dim_q
!
                  pq = dim_p*(q-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_qrsp(q,r,s,p) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_2341
!
!
   subroutine sort_1234_to_2314(x_pqrs, x_qrps, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_qrps (i.e., 1234 to 2314).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_q, dim_r, dim_p, dim_s), intent(out) :: x_qrps
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do p = 1, dim_p
            do r = 1, dim_r
               do q = 1, dim_q
!
                  x_qrps(q,r,p,s) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2314
!
!
   subroutine sort_1234_to_2134(x_pqrs, x_qprs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2134
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!
!!    Reorders the array x_pqrs to x_qprs (i.e., 1234 to 2134).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_q, dim_p, dim_r, dim_s), intent(out) :: x_qprs
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do p = 1, dim_p
               do q = 1, dim_q
!
                  x_qprs(q, p, r, s) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2134
!
!
   subroutine sort_1234_to_2143(x_pqrs, x_qpsr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2143
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!
!!    Reorders the array x_pqrs to x_qpsr (i.e., 1234 to 2143).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_q, dim_p, dim_s, dim_r), intent(out) :: x_qpsr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do r = 1, dim_r
         do s = 1, dim_s
            do p = 1, dim_p
               do q = 1, dim_q
!
                  x_qpsr(q, p, s, r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2143
!
!
   subroutine sort_1234_to_2413(x_pqrs, x_qspr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2413
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!
!!    Reorders the array x_pqrs to x_qspr (i.e., 1234 to 2413).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_q, dim_s, dim_p, dim_r), intent(out) :: x_qspr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do r = 1, dim_r
         do p = 1, dim_p
            do s = 1, dim_s
               do q = 1, dim_q
!
                  x_qspr(q, s, p, r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2413
!
!
   subroutine sort_1234_to_1324(x_pqrs, x_prqs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1324
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_pr_qs (i.e., 1234 to 1324).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_p, dim_r, dim_q, dim_s), intent(out) :: x_prqs
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do q = 1, dim_q
            do r = 1, dim_r
               do p = 1, dim_p
!
                  x_prqs(p, r, q, s) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_1324
!
!
   subroutine squareup_and_sort_1234_to_1324(x_pqrs, x_prqs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 1324
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_pr_qs (i.e., 1234 to 1324).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*(dim_r*dim_s)/2)), intent(in) :: x_pqrs
      real(dp), dimension(dim_p, dim_r, dim_q, dim_s), intent(out)       :: x_prqs
!
      integer :: p, q, r, s, pq, rs, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do s = 1, dim_s
         do q = 1, dim_q
            do r = 1, dim_r
!
               rs = dim_r*(s-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_prqs(p,r,q,s) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_1324
!
!
   subroutine sort_1234_to_2341(x_pqrs, x_qrsp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2341
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_qrsp (i.e., 1234 to 2341).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_q, dim_r, dim_s, dim_p), intent(out) :: x_qrsp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  x_qrsp(q,r,s,p) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2341
!
!
   subroutine sort_1234_to_1342(x_pqrs, x_prsq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1342
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_prsq (i.e., 1234 to 1342).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_p, dim_r, dim_s, dim_q), intent(out) :: x_prsq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do s = 1, dim_s
            do r = 1, dim_r
               do p = 1, dim_p
!
                  x_prsq(p,r,s,q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_1342
!
!
   subroutine squareup_and_sort_1234_to_1342(x_pqrs, x_prsq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 1342
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_p_r_s_q (i.e., 1234 to 1342).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*dim_r*dim_s)/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_p, dim_r, dim_s, dim_q), intent(out)     :: x_prsq
!
      integer :: p, q, r, s, rs, pq, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do q = 1, dim_q
         do s = 1, dim_s
            do r = 1, dim_r
!
               rs = dim_r*(s-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
!
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_prsq(p,r,s,q) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_1342
!
!
   subroutine sort_1234_to_1432(x_pqrs, x_psrq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1432
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_psrq (i.e., 1234 to 1432).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_p, dim_s, dim_r, dim_q), intent(out) :: x_psrq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do r = 1, dim_r
            do s = 1, dim_s
               do p = 1, dim_p
!
                  x_psrq(p,s,r,q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_1432
!
!
   subroutine squareup_and_sort_1234_to_4312(x_pqrs, x_srpq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Square up and sort 1234 to 4312
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre
!!    and Andreas Skeidsvoll, 2018
!!
!!    Squares up and reorders the array x_pqrs to x_sr_pq (i.e., 1234 to 4312).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*(dim_r*dim_s))/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_s, dim_r, dim_p, dim_q), intent(out)       :: x_srpq
!
      integer :: p, q, r, s, pq, rs, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do q = 1, dim_q
         do p = 1, dim_p
!
            pq = dim_p*(q-1) + p
!
            do r = 1, dim_r
               do s = 1, dim_s
   !
                  rs = dim_r*(s-1) + r
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
   !
                  x_srpq(s,r,p,q) = x_pqrs(pqrs)
   !
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_4312
!
!
   subroutine sort_1234_to_4312(x_pqrs, x_srpq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4312
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_srpq (i.e., 1234 to 4312).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_s, dim_r, dim_p, dim_q), intent(out) :: x_srpq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(p,q,r,s)
      do q = 1, dim_q
         do p = 1, dim_p
            do r = 1, dim_r
               do s = 1, dim_s
!
                  x_srpq(s,r,p,q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_4312
!
!
   subroutine sort_1234_to_1423(x_pqrs, x_psqr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1423
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_psqr (i.e., 1234 to 1423).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_p, dim_s, dim_q, dim_r), intent(out) :: x_psqr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do r = 1, dim_r
         do q = 1, dim_q
            do s = 1, dim_s
               do p = 1, dim_p
!
                  x_psqr(p,s,q,r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_1423
!
!
   subroutine add_1423_to_1234(scalar, x_psqr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1423 to 1234
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and Andreas Skeidsvoll, 2018
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_psqr(p,s,q,r)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_p, dim_s, dim_q, dim_r), intent(in)    :: x_psqr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_psqr(p,s,q,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_1423_to_1234
!
!
   subroutine add_1432_to_1234(scalar, x_psrq, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1432 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x(p,s,r,q)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_p, dim_s, dim_r, dim_q), intent(in)    :: x_psrq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_psrq(p,s,r,q)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_1432_to_1234
!
!
   subroutine add_1342_to_1234(scalar, x_prsq, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1342 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_prsq(p,r,s,q)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_p, dim_r, dim_s, dim_q), intent(in)    :: x_prsq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_prsq(p,r,s,q)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_1342_to_1234
!
!
   subroutine add_1324_to_1234(scalar, x_prqs, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1342 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Andreas Skeidsvoll, Apr 2019: Created by modifying add_1342_to_1234
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_prqs(p,r,q,s)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_p, dim_r, dim_q, dim_s), intent(in)    :: x_prqs
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_prqs(p,r,q,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_1324_to_1234
!
!
   subroutine add_1243_to_1234(scalar, x_pqsr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1243 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_pqsr(p,q,s,r)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_p, dim_q, dim_s, dim_r), intent(in)    :: x_pqsr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_pqsr(p,q,s,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_1243_to_1234
!
!
   subroutine add_3412_to_1234(scalar, x_rspq, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 3412 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x(r,s,p,q)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_r, dim_s, dim_p, dim_q), intent(in)    :: x_rspq
!
      call add_21_to_12(scalar, x_rspq, y_pqrs, dim_p*dim_q, dim_r*dim_s)
!
   end subroutine add_3412_to_1234
!
!
   subroutine add_3421_to_1234(scalar, x_rsqp, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 3421 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_rsqp(r,s,q,p)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_r, dim_s, dim_q, dim_p), intent(in)    :: x_rsqp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(p,q,r,s)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_rsqp(r,s,q,p)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_3421_to_1234
!
!
   subroutine add_2341_to_1234(scalar, x_qrsp, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 2341 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_qrsp(q,r,s,p)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_q, dim_r, dim_s, dim_p), intent(in)    :: x_qrsp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_qrsp(q,r,s,p)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_2341_to_1234
!
!
   subroutine add_2143_to_1234(scalar, x_qpsr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 2143 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x(q,p,s,r)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_q, dim_p, dim_s, dim_r), intent(in)    :: x_qpsr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_qpsr(q,p,s,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_2143_to_1234
!
!
   subroutine add_2134_to_1234(scalar, x_qprs, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 2143 to 1234
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_qprs(q,p,r,s)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_q, dim_p, dim_r, dim_s), intent(in)    :: x_qprs
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_qprs(q,p,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_2134_to_1234
!
!
   subroutine add_3214_to_1234(scalar, x_rqps, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 3214 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x(r,q,p,s)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_r, dim_q, dim_p, dim_s), intent(in)    :: x_rqps
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_rqps(r,q,p,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_3214_to_1234
!
!
   subroutine add_4231_to_1234(scalar, x_sqrp, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4231 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x(s,q,r,p)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_s, dim_q, dim_r, dim_p), intent(in)    :: x_sqrp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_sqrp(s,q,r,p)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_4231_to_1234
!
!
subroutine add_2413_to_1234(scalar, x_qspr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 2431 to 1234
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Alexander Paul, Jan 2019
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_qspr(q,s,p,r)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_q, dim_s, dim_p, dim_r), intent(in)    :: x_qspr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_qspr(q,s,p,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_2413_to_1234
!
!
   subroutine add_2431_to_1234(scalar, x_qsrp, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 2431 to 1234
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_qsrp(q,s,r,p)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_q, dim_s, dim_r, dim_p), intent(in)    :: x_qsrp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_qsrp(q,s,r,p)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_2431_to_1234
!
!
   subroutine add_4213_to_1234(scalar, x_sqpr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4213 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_sqpr(s,q,p,r)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_s, dim_q, dim_p, dim_r), intent(in)    :: x_sqpr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_sqpr(s,q,p,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_4213_to_1234
!
!
   subroutine add_4321_to_1234(scalar, x_srqp, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4321 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_srqp(s,r,q,p)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_s, dim_r, dim_q, dim_p), intent(in)    :: x_srqp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_srqp(s,r,q,p)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_4321_to_1234
!
!
   subroutine sort_1234_to_4321(x_pqrs, x_srqp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4321
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and Andreas Skeidsvoll, 2018
!!
!!    Reorders the array x_pqrs to x_srqp (i.e., 1234 to 4321).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_s, dim_r, dim_q, dim_p), intent(out) :: x_srqp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do p = 1, dim_p
         do q = 1, dim_q
            do r = 1, dim_r
               do s = 1, dim_s
   !
                  x_srqp(s,r,q,p) = x_pqrs(p,q,r,s)
   !
               enddo
            enddo
         enddo
      enddo
   !$omp end parallel do
!
   end subroutine sort_1234_to_4321
!
!
   subroutine add_4312_to_1234(scalar, x_srpq, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4312 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_srpq(s,r,p,q)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_s, dim_r, dim_p, dim_q), intent(in)    :: x_srpq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_srpq(s,r,p,q)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_4312_to_1234
!
!
   subroutine add_4123_to_1234(scalar, x_spqr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4123 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_spqr(s,p,q,r)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_s, dim_p, dim_q, dim_r), intent(in)    :: x_spqr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_spqr(s,p,q,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_4123_to_1234
!
!
   subroutine add_4132_to_1234(scalar, x_sprq, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4132 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_sprq(s,p,r,q)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_s, dim_p, dim_r, dim_q), intent(in)    :: x_sprq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_sprq(s,p,r,q)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_4132_to_1234
!
!
   subroutine squareup_and_sort_1234_to_1432(x_pqrs, x_psrq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 1432
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the packed array x_pqrs to unpacked x_psrq (i.e., 1234 to 1432).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*(dim_r*dim_s))/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_p, dim_s, dim_r, dim_q), intent(out)       :: x_psrq
!
      integer :: p, q, r, s, pq, rs, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do q = 1, dim_q
         do r = 1, dim_r
            do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
!
               do s = 1, dim_s
!
                  rs = dim_r*(s-1) + r
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_psrq(p,s,r,q) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_1432
!
!
   subroutine squareup_and_sort_1234_to_1243(x_pqrs, x_pqsr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 1243
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the packed array x_pqrs to unpacked x_pqsr (i.e., 1234 to 1432).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*(dim_r*dim_s))/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_p, dim_q, dim_s, dim_r), intent(out)       :: x_pqsr
!
      integer :: p, q, r, s, pq, rs, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do r = 1, dim_r
         do s = 1, dim_s
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_pqsr(p,q,s,r) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_1243
!
!
   subroutine squareup_and_sort_1234_to_1423(x_pqrs, x_psqr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 1423
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the packed array x_pqrs to unpacked x_psqr (i.e., 1234 to 1423).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*(dim_r*dim_s))/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_p, dim_s, dim_q, dim_r), intent(out)       :: x_psqr
!
      integer :: p, q, r, s, pq, rs, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do r = 1, dim_r
         do q = 1, dim_q
            do s = 1, dim_s
!
               rs = dim_r*(s-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_psqr(p,s,q,r) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_1423
!
!
   subroutine sort_1234_to_3214(x_pqrs, x_rqps, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3214
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_rq_ps (i.e., 1234 to 3214).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_r, dim_q, dim_p, dim_s), intent(out) :: x_rqps
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do p = 1, dim_p
            do q = 1, dim_q
               do r = 1, dim_r
!
                  x_rqps(r,q,p,s) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_3214
!
!
   subroutine sort_1234_to_4231(x_pqrs, x_sqrp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4231
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_sqrp (i.e., 1234 to 4231).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_s, dim_q, dim_r, dim_p), intent(out) :: x_sqrp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do p = 1, dim_p
         do r = 1, dim_r
            do q = 1, dim_q
               do s = 1, dim_s
!
                  x_sqrp(s,q,r,p) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_4231
!
!
   subroutine sort_1234_to_4213(x_pqrs, x_sqpr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4213
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_sq_pr (i.e., 1234 to 4213).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_s, dim_q, dim_p, dim_r), intent(out) :: x_sqpr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do r = 1, dim_r
         do p = 1, dim_p
            do q = 1, dim_q
               do s = 1, dim_s
!
                  x_sqpr(s,q,p,r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_4213
!
!
   subroutine squareup_and_sort_1234_to_4213(x_pqrs, x_sqpr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Square up and sort 1234 to 4213
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre
!!    and Andreas Skeidsvoll, 2018
!!
!!    Reorders the array x_pqrs to x_sqpr (i.e., 1234 to 4213).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*dim_r*dim_s)/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_s, dim_q, dim_p, dim_r), intent(out)     :: x_sqpr
!
      integer :: p, q, r, s, pq, rs, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do r = 1, dim_r
         do p = 1, dim_p
            do q = 1, dim_q
!
               pq = dim_p*(q-1) + p
!
               do s = 1, dim_s
!
                  rs = dim_r*(s-1) + r
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_sqpr(s,q,p,r) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_4213
!
!
   subroutine squareup_and_sort_1234_to_3214(x_pqrs, x_rqps, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 3214
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the packed array x_pqrs to unpacked x_rqps (i.e., 1234 to 3214).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*dim_r*dim_s)/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_r, dim_q, dim_p, dim_s), intent(out)     :: x_rqps
!
      integer :: p, q, r, s, pq, rs, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do s = 1, dim_s
         do p = 1, dim_p
            do q = 1, dim_q
!
               pq = dim_p*(q-1) + p
!
               do r = 1, dim_r
!
                  rs = dim_r*(s-1) + r
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_rqps(r,q,p,s) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_3214
!
!
   subroutine sort_1234_to_3241(x_pqrs, x_rqsp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3241
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_rqsp (i.e., 1234 to 3241).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_r, dim_q, dim_s, dim_p), intent(out) :: x_rqsp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do p = 1, dim_p
         do s = 1, dim_s
            do q = 1, dim_q
               do r = 1, dim_r
!
                  x_rqsp(r,q,s,p) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_3241
!
!
   subroutine sort_1234_to_1243(x_pqrs, x_pqsr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1243
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_pqsr (i.e., 1234 to 1243).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_p, dim_q, dim_s, dim_r), intent(out) :: x_pqsr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(p, q, r, s)
      do r = 1,dim_r
         do s = 1,dim_s
            do q = 1, dim_q
               do p = 1, dim_p
!
                  x_pqsr(p,q,s,r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_1243
!
!
!     -::- Squareup and packin and related routines -::-
!     --------------------------------------------------
!
   subroutine add_to_packed(packed, unpacked, N)
!!
!!    Add to packed
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Adds a symmetric unpacked N x N matrix to a packed N x N
!!    matrix.
!!
      implicit none
!
      integer, intent(in) :: N
!
      real(dp), dimension(N*(N+1)/2), intent(inout) :: packed
      real(dp), dimension(N, N), intent(in)         :: unpacked
!
      integer :: i, j, ij
!
!$omp parallel do schedule(static) private(i,j, ij)
      do i = 1, N
         do j = i, N
            ij = (max(i,j)*(max(i,j)-3)/2) + i + j
            packed(ij) = packed(ij) + unpacked(i,j)
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_to_packed
!
!
   subroutine symmetrize_and_add_to_packed(packed, unpacked, N)
!!
!!    Symmetrize and add to packed
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Symmetrizes unpacked N x N matrix and adds the symmtrized sum
!!    to a packed N x N matrix.
!!
      implicit none
!
      integer, intent(in) :: N
!
      real(dp), dimension(N*(N+1)/2), intent(out) :: packed
      real(dp), dimension(N,N), intent(in)        :: unpacked
!
      integer :: i, j, ij
!
!$omp parallel do schedule(static) private(i,j, ij)
      do i = 1, N
         do j = i, N
            ij = (max(i,j)*(max(i,j)-3)/2) + i + j
            packed(ij) = packed(ij) + unpacked(i,j) + unpacked(j,i)
         enddo
      enddo
!$omp end parallel do
!
   end subroutine symmetrize_and_add_to_packed
!
!
   subroutine squareup_anti(packed, unpacked, N)
!!
!!    Square up packed antisymmetric matrix
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Squares up to full dimension (N x N) of packed matrix. The packed
!!    antisymmetric matrix contains the strictly lower triangular part
!     of the full unpacked matrix.
!!
      implicit none
!
      integer, intent(in) :: N
!
      real(dp), dimension(N*(N-1)/2), intent(in)  :: packed
      real(dp), dimension(N,N), intent(out)       :: unpacked
!
      integer :: i, j
!
!     Set diagonal to zero
!
      do i = 1, N
!
         unpacked(i, i) = zero
!
      enddo
!
!     Set lower and upper strictly triangular parts
!
      do i = 2, N
         do j = 1, i - 1
!
            unpacked(i, j) = packed((max(i-1,j)*(max(i-1,j)-3)/2) + i-1 + j)
            unpacked(j, i) = -packed((max(i-1,j)*(max(i-1,j)-3)/2) + i-1 + j)
!
         enddo
      enddo
!
   end subroutine squareup_anti
!
!
   subroutine squareup(packed,unpacked,N)
!!
!!    Square up packed symmetric matrix
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Squares up to full dimension (N x N) of packed matrices.
!!
      implicit none
!
      integer, intent(in) :: N
!
      real(dp), dimension(N*(N+1)/2), intent(in) :: packed
      real(dp), dimension(N,N), intent(out)      :: unpacked
!
      integer :: i, j
!
!$omp parallel do schedule(static) private(i,j)
      do j = 1, N
         do i = 1, N
            unpacked(i, j) = packed((max(i,j)*(max(i,j)-3)/2) + i + j)
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup
!
!
   subroutine packin(packed,unpacked,N)
!!
!!    Pack in symmetric matrix
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Pack down full square matrix of dimension N x N.
!!
      implicit none
!
      integer, intent(in) :: N
!
      real(dp), dimension(N*(N+1)/2), intent(out) :: packed
      real(dp), dimension(n,N), intent(in)        :: unpacked
!
      integer :: i, j
!
!$omp parallel do schedule(static) private(i,j)
      do i = 1, N
         do j = 1, i
!
            packed((max(i,j)*(max(i,j)-3)/2) + i + j) = unpacked(i, j)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine packin
!
!
   subroutine packin_anti(packed, unpacked, N)
!!
!!    Pack in anti-symmetric matrix
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Pack down full square anti-symmetric matrix of dimension N x N,
!!    where the strictly lower triangular part of the unpacked matrix
!!    is stored in packed form.
!!
      implicit none
!
      integer, intent(in) :: N
!
      real(dp), dimension(N*(N-1)/2), intent(out) :: packed
      real(dp), dimension(N,N), intent(in)        :: unpacked
!
      integer :: i, j
!
      do i = 2, N
         do j = 1, i - 1
!
            packed((max(i-1,j)*(max(i-1,j)-3)/2) + i-1 + j) = unpacked(i, j)
!
         enddo
      enddo
!
   end subroutine packin_anti
!
!
end module reordering
