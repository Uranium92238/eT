!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module packed_array_utilities_r
!
!!
!!    Real packed array utilities module
!!
!!    Routines that perform various operations on packed and double packed arrays.
!!
!!    Some definitions:
!!
!!    Packed array: 
!!    Symmetric square nxn matrix M
!!    where only the upper triangular matrix is stored 
!!    in an array P of size n*(n+1)/2. 
!!    P(p + q*(q-1)) = M(p,q) with p >= q
!!
!!    Double packed:
!!    4d arrays where where the indices are combined into packed compound indices.
!!    If the 4d array has dimensions M(dim1, dim1, dim2, dim2), the double packed array
!!    has dimensions P(dim1*(dim1+1)/2, dim2*(dim2+1)/2)
!!    
!!    Single packed:
!!    Similar to above, but only a single pair of indices are compounded and packed.
!!    M(dim1, dim1, dim2, dim3) -> P(dim1*(dim1+1)/2, dim2, dim3)
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
!!
!
      use kinds
!
contains
!
!
   subroutine scale_double_packed_diagonal(x, dim_p, dim_q, factor)
!!
!!    Scale double packed diagonal
!!    Written by Rolf H. Myhre Jan. 2021
!!
!!    x :: double packed 4d array
!!
!!    x(pp,qq) = factor*x(pp,qq)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q
!
      real(dp), dimension(dim_p*(dim_p+1)/2, dim_q*(dim_q+1)/2), intent(inout) :: x
!
      real(dp), intent(in) :: factor
!
      integer :: p, q
!
!$omp parallel do private(p,q)
      do p = 1, dim_p
         do q = 1, dim_q
            x(p*(p+1)/2, q*(q+1)/2) = factor*x(p*(p+1)/2, q*(q+1)/2)
         enddo
      enddo
!$omp end parallel do
!
   end subroutine scale_double_packed_diagonal
!
!
   subroutine scale_double_packed_blocks(x, dim_p, dim_q, q_factor, p_factor)
!!
!!    Scale double packed blocks
!!    Written by Rolf H. Myhre Jan. 2021
!!
!!    Scale diagonal blocks of double packed array x
!!
!!    x(:,qq) => q_factor*x(:,qq)
!!    x(pp,:) => p_factor*x(pp,:)
!!
!!    p_factor and q_factor are optional and will be skipped if not present
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q
!
      real(dp), dimension(dim_p*(dim_p+1)/2, dim_q*(dim_q+1)/2), intent(inout) :: x
!
      real(dp), optional, intent(in) :: q_factor, p_factor
!
      integer :: p, q
!
      if(present(q_factor)) then
!$omp parallel do private(q)
         do q = 1, dim_q
            x(:, q*(q+1)/2) = q_factor*x(:, q*(q+1)/2)
         enddo
!$omp end parallel do
      endif
!
      if(present(p_factor)) then
!$omp parallel do private(p)
         do p = 1, dim_p
            x(p*(p+1)/2, :) = p_factor*x(p*(p+1)/2, :)
         enddo
!$omp end parallel do
      endif
!
   end subroutine scale_double_packed_blocks
!
!
   subroutine construct_plus_minus_2413_from_packed(x, x_p, x_m, dim_1, dim_2)
!!
!!    Construct plus minus 2413 from packed
!!    Written by Rolf H. Myhre, Eirik F. Kjønstad, and Sarai D. Folkestad Jan. 2021
!!
!!    x   :: packed 4d array
!!    x_p :: double packed plus combination,  x_p(qs,pr) = x(p,q,r,s) + x(p,s,r,q)
!!    x_m :: double packed minus combination, x_m(qs,pr) = x(p,q,r,s) - x(p,s,r,q)
!!
      implicit none
!
      integer, intent(in) :: dim_1, dim_2
!
      real(dp), dimension((dim_1*dim_2)*(dim_1*dim_2+1)/2), intent(in) :: x
!
      real(dp), dimension(dim_1*(dim_1+1)/2, dim_2*(dim_2+1)/2), intent(out) :: x_p
      real(dp), dimension(dim_1*(dim_1-1)/2, dim_2*(dim_2-1)/2), intent(out) :: x_m
!
      integer :: p, r, q, s, pr_p, qs_p, pr_m, qs_m
      integer :: pq, rs, ps, rq, pqrs, psrq
!
!$omp parallel do private(p,q,r,s,pq,rs,rq,ps,pqrs,psrq,pr_p,qs_p,pr_m,qs_m) &
!$omp schedule(guided)
      do r = 1, dim_2
         do p = 1, r
!
            pr_p = r*(r-1)/2 + p
            pr_m = (r-1)*(r-2)/2 + p
!
            do s = 1, dim_1
!
               ps = dim_2*(s-1) + p
               rs = dim_2*(s-1) + r
!
               do q = 1, s
!
                  qs_p = s*(s-1)/2 + q
                  qs_m = (s-1)*(s-2)/2 + q
!
                  pq = dim_2*(q-1) + p
                  rq = dim_2*(q-1) + r
!
                  pqrs = rs*(rs-1)/2 + pq
                  psrq = max(rq,ps)*(max(rq,ps)-3)/2 + rq + ps
!
                  x_p(qs_p,pr_p) = x(pqrs) + x(psrq)
                  if(p .ne. r .and. q .ne. s) then
                     x_m(qs_m,pr_m) = x(pqrs) - x(psrq)
                  endif
!
               enddo
            enddo
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_plus_minus_2413_from_packed
!
!
   subroutine construct_plus_minus_1324_from_RFP(x, x_p, x_m, dim_1, dim_2)
!!
!!    Construct plus minus 1324 from RFP
!!    Written by Rolf H. Myhre, Eirik F. Kjønstad, and Sarai D. Folkestad Jan. 2021
!!
!!    x   :: packed 4d array in the full rectangular packed N U format
!!    x_p :: double packed plus combination,  x_p(pr,qs) = x(p,q,r,s) + x(p,s,r,q)
!!    x_m :: double packed minus combination, x_m(pr,qs) = x(p,q,r,s) - x(p,s,r,q)
!!
      implicit none
!
      integer, intent(in) :: dim_1, dim_2
!
      real(dp), dimension(dim_1*dim_2 + mod(dim_1*dim_2+1,2), (dim_1*dim_2+1)/2), intent(in) :: x
!
      real(dp), dimension(dim_1*(dim_1+1)/2, dim_2*(dim_2+1)/2), intent(out) :: x_p
      real(dp), dimension(dim_1*(dim_1-1)/2, dim_2*(dim_2-1)/2), intent(out) :: x_m
!
      integer :: tridim
      integer :: p,q,r,s
      integer :: pr_p, pr_m, qs_p, qs_m
      integer :: pq, rs, ps, rq
      integer :: pq_rfp, rs_rfp, ps_rfp, rq_rfp
!
      tridim  = dim_1*dim_2/2
!
!$omp parallel do private(p,q,r,s,pq,rs,ps,rq,pq_rfp,rs_rfp,ps_rfp,rq_rfp,pr_p,qs_p,pr_m,qs_m) &
!$omp schedule(guided)
      do s = 1, dim_2
         do q = 1, s
!
            qs_p = q + s*(s-1)/2
            qs_m = q + (s-1)*(s-2)/2
!
            do r = 1, dim_1
!
               rq = r + (q-1)*dim_1
               rs = r + (s-1)*dim_1
!
               do p = 1, r
!
                  pr_p = p + r*(r-1)/2
                  pr_m = p + (r-1)*(r-2)/2
!
                  pq = p + (q-1)*dim_1
                  ps = p + (s-1)*dim_1
!
                  if (rs .gt. tridim) then
                     pq_rfp = pq
                     rs_rfp = rs - tridim
                  else
                     pq_rfp = rs + tridim + 1
                     rs_rfp = pq
                  endif
!
                  if (ps .ge. rq) then
                     if (ps .gt. tridim) then
                        ps_rfp = rq
                        rq_rfp = ps - tridim
                     else
                        ps_rfp = ps + tridim + 1
                        rq_rfp = rq
                     endif
                  else
                     if (rq .gt. tridim) then
                        ps_rfp = ps
                        rq_rfp = rq - tridim
                     else
                        ps_rfp = rq + tridim + 1
                        rq_rfp = ps
                     endif
                  endif
!
                  x_p(pr_p, qs_p) = x(pq_rfp, rs_rfp) + x(ps_rfp, rq_rfp)
!
                  if (p .ne. r .and. q .ne. s) then
!
                     x_m(pr_m, qs_m) = x(pq_rfp, rs_rfp) - x(ps_rfp, rq_rfp)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_plus_minus_1324_from_RFP
!
!
   subroutine construct_plus_minus_1324_from_full(x, x_p, x_m, dim_1, dim_2, dim_3)
!!
!!    Construct plus minus 1324 from full
!!    Written by Rolf H. Myhre, Eirik F. Kjønstad, and Sarai D. Folkestad Jan. 2021
!!
!!    x   :: full 4d array
!!    x_p :: single packed plus combination,  x_p(pr,q,s) = x(p,q,r,s) + x(p,s,r,q)
!!    x_m :: single packed minus combination, x_m(pr,q,s) = x(p,q,r,s) - x(p,s,r,q)
!!
      implicit none
!
      integer, intent(in) :: dim_1, dim_2, dim_3
!
      real(dp), dimension(dim_1, dim_2, dim_1, dim_3), intent(in) :: x
!
      real(dp), dimension(dim_1*(dim_1+1)/2, dim_2, dim_3), intent(out) :: x_p
      real(dp), dimension(dim_1*(dim_1-1)/2, dim_2, dim_3), intent(out) :: x_m
!
      integer :: p, r, q, s, pr_p, pr_m
!
!$omp parallel do private(p,q,r,s,pr_p,pr_m), collapse(2)
      do s = 1, dim_3
         do q = 1, dim_2
            do r = 1, dim_1
               do p = 1, r
!
                  pr_p = p + r*(r-1)/2
                  pr_m = p + (r-1)*(r-2)/2
!
                  x_p(pr_p,q,s) = x(p,q,r,s) + x(r,q,p,s)
!
                  if(p .ne. r) then
                     x_m(pr_m,q,s) = x(p,q,r,s) - x(r,q,p,s)
                  endif
!
               enddo
            enddo
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_plus_minus_1324_from_full
!
!
   subroutine add_double_packed_plus_minus(x, x_p, x_m, dim_1, dim_2, dim_full, offset)
!!
!!    Add double packed plus minus
!!    Written by Rolf H. Myhre, Eirik F. Kjønstad, and Sarai D. Folkestad Jan. 2021
!!
!!    x   :: packed 4d array
!!    x_p :: double packed plus combination 
!!    x_m :: double packed minus combination
!!
!!    x(pqrs) = x(pqrs) + x_p(qs,pr) + x_m(qs,pr)
!!    x(psrq) = x(psrq) + x_p(qs,pr) - x_m(qs,pr)
!!
      implicit none
!
      integer, intent(in) :: dim_1, dim_2, dim_full, offset
!
      real(dp), dimension((dim_full*dim_1)*(dim_full*dim_1+1)/2), intent(inout) :: x
!
      real(dp), dimension(dim_1*(dim_1+1)/2, dim_2*(dim_2+1)/2), intent(in) :: x_p
      real(dp), dimension(dim_1*(dim_1-1)/2, dim_2*(dim_2-1)/2), intent(in) :: x_m
!
      integer :: p, r, q, s
      integer :: pr_p, pr_m, qs_p, qs_m
      integer :: pq, rs, ps, rq, pqrs, psrq
!
!$omp parallel do private(p,q,r,s,pq,rs,rq,ps,pqrs,psrq,pr_p,qs_p,pr_m,qs_m) &
!$omp schedule(guided)
      do s = 1, dim_1
         do q = 1, s
!
            qs_p = q + s*(s-1)/2
            qs_m = q + (s-1)*(s-2)/2
!
            do r = 1, dim_2
!
               rs = r + dim_full*(s-1) + offset
               rq = r + dim_full*(q-1) + offset
!
               do p = 1, r
!
                  pr_p = p + r*(r-1)/2
                  pr_m = p + (r-1)*(r-2)/2
!
                  pq = p + dim_full*(q-1) + offset
                  ps = p + dim_full*(s-1) + offset
!
                  pqrs = pq + rs*(rs-1)/2
                  psrq = max(rq,ps)*(max(rq,ps)-3)/2 + rq + ps
!
                  if (q .ne. s .and. p .ne. r) then
                     x(pqrs) = x(pqrs) + x_p(qs_p, pr_p) + x_m(qs_m, pr_m)
                     x(psrq) = x(psrq) + x_p(qs_p, pr_p) - x_m(qs_m, pr_m)
                  else
                     x(pqrs) = x(pqrs) + x_p(qs_p, pr_p)
                     x(psrq) = x(psrq) + x_p(qs_p, pr_p)
                  endif
!
               enddo
            enddo
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_double_packed_plus_minus
!
!
   subroutine add_single_packed_plus_minus(x, x_p, x_m, &
                                           dim_1, dim_2, dim_3, &
                                           dim_full, offset_2, offset_3)
!!
!!    Add single packed plus minus
!!    Written by Rolf H. Myhre, Eirik F. Kjønstad, and Sarai D. Folkestad Jan. 2021
!!
!!    x   :: packed 4d array
!!    x_p :: single packed plus combination
!!    x_m :: single packed minus combination
!!
!!    x(pqrs) = x(pqrs) + x_p(qs,p,r) + x_m(qs,p,r)
!!    x(psrq) = x(psrq) + x_p(qs,p,r) - x_m(qs,p,r)
!!
      implicit none
!
      integer, intent(in) :: dim_1, dim_2, dim_3, dim_full, offset_2, offset_3
!
      real(dp), dimension((dim_full*dim_1)*(dim_full*dim_1+1)/2), intent(inout) :: x
!
      real(dp), dimension(dim_1*(dim_1+1)/2, dim_2, dim_3), intent(in) :: x_p
      real(dp), dimension(dim_1*(dim_1-1)/2, dim_2, dim_3), intent(in) :: x_m
!
      integer :: p, r, q, s, qs_p, qs_m
      integer :: pq, rs, pqrs, ps, rq, psrq
!
!$omp parallel do private(p,q,r,s,pq,rs,pqrs,ps,rq,psrq,qs_p,qs_m) &
!$omp schedule(guided)
      do s = 1, dim_1
         do q = 1, s
!
            qs_p = q + s*(s-1)/2
            qs_m = q + (s-1)*(s-2)/2
!
            do r = 1, dim_3
!
               rq = r + dim_full*(q-1) + offset_3
               rs = r + dim_full*(s-1) + offset_3
!
               do p = 1, dim_2
!
                  pq = p + dim_full*(q-1) + offset_2
                  ps = p + dim_full*(s-1) + offset_2
!
                  pqrs = max(pq,rs)*(max(pq,rs)-3)/2 + pq + rs
                  psrq = max(ps,rq)*(max(ps,rq)-3)/2 + ps + rq
!
                  if(q .ne. s) then
                     x(pqrs) = x(pqrs) + x_p(qs_p, p, r) + x_m(qs_m, p, r)
                     x(psrq) = x(psrq) + x_p(qs_p, p, r) - x_m(qs_m, p, r)
                  else
                     x(pqrs) = x(pqrs) + x_p(qs_p, p, r)
                     x(psrq) = x(psrq) + x_p(qs_p, p, r)
                  endif
!
               enddo
            enddo
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_single_packed_plus_minus
!
!
end module packed_array_utilities_r
!
!
module packed_array_utilities_c
!
!!
!!    Complex packed array utilities module
!!
!!    Routines that perform various operations on packed and double packed arrays.
!!
!!    Some definitions:
!!
!!    Packed array: 
!!    Symmetric scuare nxn matrix M
!!    where only the upper triangular matrix is stored 
!!    in an array P of size n*(n+1)/2. 
!!    P(p + q*(q-1)) = M(p,q) with p >= q
!!
!!    Double packed:
!!    4d arrays where where the indices are combined into packed compound indices.
!!    If the 4d array has dimensions M(dim1, dim1, dim2, dim2), the double packed array
!!    has dimensions P(dim1*(dim1+1)/2, dim2*(dim2+1)/2)
!!    
!!    Single packed:
!!    Similar to above, but only a single pair of indices are compounded and packed.
!!    M(dim1, dim1, dim2, dim3) -> P(dim1*(dim1+1)/2, dim2, dim3)
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
      use kinds
!
contains
!
!
   subroutine scale_double_packed_diagonal(x, dim_p, dim_q, factor)
!!
!!    Scale double packed diagonal
!!    Written by Rolf H. Myhre Jan. 2021
!!
!!    x :: double packed 4d array
!!
!!    x(pp,qq) = factor*x(pp,qq)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q
!
      complex(dp), dimension(dim_p*(dim_p+1)/2, dim_q*(dim_q+1)/2), intent(inout) :: x
!
      complex(dp), intent(in) :: factor
!
      integer :: p, q
!
!$omp parallel do private(p,q)
      do p = 1, dim_p
         do q = 1, dim_q
            x(p*(p+1)/2, q*(q+1)/2) = factor*x(p*(p+1)/2, q*(q+1)/2)
         enddo
      enddo
!$omp end parallel do
!
   end subroutine scale_double_packed_diagonal
!
!
   subroutine scale_double_packed_blocks(x, dim_p, dim_q, q_factor, p_factor)
!!
!!    Scale double packed blocks
!!    Written by Rolf H. Myhre Jan. 2021
!!
!!    Scale diagonal blocks of double packed array x
!!
!!    x(:,qq) => q_factor*x(:,qq)
!!    x(pp,:) => p_factor*x(pp,:)
!!
!!    p_factor and q_factor are optional and will be skipped if not present
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q
!
      complex(dp), dimension(dim_p*(dim_p+1)/2, dim_q*(dim_q+1)/2), intent(inout) :: x
!
      complex(dp), optional, intent(in) :: q_factor, p_factor
!
      integer :: p, q
!
      if(present(q_factor)) then
!$omp parallel do private(q)
         do q = 1, dim_q
            x(:, q*(q+1)/2) = q_factor*x(:, q*(q+1)/2)
         enddo
!$omp end parallel do
      endif
!
      if(present(p_factor)) then
!$omp parallel do private(p)
         do p = 1, dim_p
            x(p*(p+1)/2, :) = p_factor*x(p*(p+1)/2, :)
         enddo
!$omp end parallel do
      endif
!
   end subroutine scale_double_packed_blocks
!
!
   subroutine construct_plus_minus_2413_from_packed(x, x_p, x_m, dim_1, dim_2)
!!
!!    Construct plus minus 2413 from packed
!!    Written by Rolf H. Myhre, Eirik F. Kjønstad, and Sarai D. Folkestad Jan. 2021
!!
!!    x   :: packed 4d array
!!    x_p :: double packed plus combination,  x_p(qs,pr) = x(p,q,r,s) + x(p,s,r,q)
!!    x_m :: double packed minus combination, x_m(qs,pr) = x(p,q,r,s) - x(p,s,r,q)
!!
      implicit none
!
      integer, intent(in) :: dim_1, dim_2
!
      complex(dp), dimension((dim_1*dim_2)*(dim_1*dim_2+1)/2), intent(in) :: x
!
      complex(dp), dimension(dim_1*(dim_1+1)/2, dim_2*(dim_2+1)/2), intent(out) :: x_p
      complex(dp), dimension(dim_1*(dim_1-1)/2, dim_2*(dim_2-1)/2), intent(out) :: x_m
!
      integer :: p, r, q, s, pr_p, qs_p, pr_m, qs_m
      integer :: pq, rs, ps, rq, pqrs, psrq
!
!$omp parallel do private(p,q,r,s,pq,rs,rq,ps,pqrs,psrq,pr_p,qs_p,pr_m,qs_m) &
!$omp schedule(guided)
      do r = 1, dim_2
         do p = 1, r
!
            pr_p = r*(r-1)/2 + p
            pr_m = (r-1)*(r-2)/2 + p
!
            do s = 1, dim_1
!
               ps = dim_2*(s-1) + p
               rs = dim_2*(s-1) + r
!
               do q = 1, s
!
                  qs_p = s*(s-1)/2 + q
                  qs_m = (s-1)*(s-2)/2 + q
!
                  pq = dim_2*(q-1) + p
                  rq = dim_2*(q-1) + r
!
                  pqrs = rs*(rs-1)/2 + pq
                  psrq = max(rq,ps)*(max(rq,ps)-3)/2 + rq + ps
!
                  x_p(qs_p,pr_p) = x(pqrs) + x(psrq)
                  if(p .ne. r .and. q .ne. s) then
                     x_m(qs_m,pr_m) = x(pqrs) - x(psrq)
                  endif
!
               enddo
            enddo
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_plus_minus_2413_from_packed
!
!
   subroutine construct_plus_minus_1324_from_RFP(x, x_p, x_m, dim_1, dim_2)
!!
!!    Construct plus minus 1324 from RFP
!!    Written by Rolf H. Myhre, Eirik F. Kjønstad, and Sarai D. Folkestad Jan. 2021
!!
!!    x   :: packed 4d array in the full rectangular packed N U format
!!    x_p :: double packed plus combination,  x_p(pr,qs) = x(p,q,r,s) + x(p,s,r,q)
!!    x_m :: double packed minus combination, x_m(pr,qs) = x(p,q,r,s) - x(p,s,r,q)
!!
      implicit none
!
      integer, intent(in) :: dim_1, dim_2
!
      complex(dp), dimension(dim_1*dim_2 + mod(dim_1*dim_2+1,2), (dim_1*dim_2+1)/2), intent(in) :: x
!
      complex(dp), dimension(dim_1*(dim_1+1)/2, dim_2*(dim_2+1)/2), intent(out) :: x_p
      complex(dp), dimension(dim_1*(dim_1-1)/2, dim_2*(dim_2-1)/2), intent(out) :: x_m
!
      integer :: tridim
      integer :: p,q,r,s
      integer :: pr_p, pr_m, qs_p, qs_m
      integer :: pq, rs, ps, rq
      integer :: pq_rfp, rs_rfp, ps_rfp, rq_rfp
!
      tridim  = dim_1*dim_2/2
!
!$omp parallel do private(p,q,r,s,pq,rs,ps,rq,pq_rfp,rs_rfp,ps_rfp,rq_rfp,pr_p,qs_p,pr_m,qs_m) &
!$omp schedule(guided)
      do s = 1, dim_2
         do q = 1, s
!
            qs_p = q + s*(s-1)/2
            qs_m = q + (s-1)*(s-2)/2
!
            do r = 1, dim_1
!
               rq = r + (q-1)*dim_1
               rs = r + (s-1)*dim_1
!
               do p = 1, r
!
                  pr_p = p + r*(r-1)/2
                  pr_m = p + (r-1)*(r-2)/2
!
                  pq = p + (q-1)*dim_1
                  ps = p + (s-1)*dim_1
!
                  if (rs .gt. tridim) then
                     pq_rfp = pq
                     rs_rfp = rs - tridim
                  else
                     pq_rfp = rs + tridim + 1
                     rs_rfp = pq
                  endif
!
                  if (ps .ge. rq) then
                     if (ps .gt. tridim) then
                        ps_rfp = rq
                        rq_rfp = ps - tridim
                     else
                        ps_rfp = ps + tridim + 1
                        rq_rfp = rq
                     endif
                  else
                     if (rq .gt. tridim) then
                        ps_rfp = ps
                        rq_rfp = rq - tridim
                     else
                        ps_rfp = rq + tridim + 1
                        rq_rfp = ps
                     endif
                  endif
!
                  x_p(pr_p, qs_p) = x(pq_rfp, rs_rfp) + x(ps_rfp, rq_rfp)
!
                  if (p .ne. r .and. q .ne. s) then
!
                     x_m(pr_m, qs_m) = x(pq_rfp, rs_rfp) - x(ps_rfp, rq_rfp)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_plus_minus_1324_from_RFP
!
!
   subroutine construct_plus_minus_1324_from_full(x, x_p, x_m, dim_1, dim_2, dim_3)
!!
!!    Construct plus minus 1324 from full
!!    Written by Rolf H. Myhre, Eirik F. Kjønstad, and Sarai D. Folkestad Jan. 2021
!!
!!    x   :: full 4d array
!!    x_p :: single packed plus combination,  x_p(pr,q,s) = x(p,q,r,s) + x(p,s,r,q)
!!    x_m :: single packed minus combination, x_m(pr,q,s) = x(p,q,r,s) - x(p,s,r,q)
!!
      implicit none
!
      integer, intent(in) :: dim_1, dim_2, dim_3
!
      complex(dp), dimension(dim_1, dim_2, dim_1, dim_3), intent(in) :: x
!
      complex(dp), dimension(dim_1*(dim_1+1)/2, dim_2, dim_3), intent(out) :: x_p
      complex(dp), dimension(dim_1*(dim_1-1)/2, dim_2, dim_3), intent(out) :: x_m
!
      integer :: p, r, q, s, pr_p, pr_m
!
!$omp parallel do private(p,q,r,s,pr_p,pr_m), collapse(2)
      do s = 1, dim_3
         do q = 1, dim_2
            do r = 1, dim_1
               do p = 1, r
!
                  pr_p = p + r*(r-1)/2
                  pr_m = p + (r-1)*(r-2)/2
!
                  x_p(pr_p,q,s) = x(p,q,r,s) + x(r,q,p,s)
!
                  if(p .ne. r) then
                     x_m(pr_m,q,s) = x(p,q,r,s) - x(r,q,p,s)
                  endif
!
               enddo
            enddo
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_plus_minus_1324_from_full
!
!
   subroutine add_double_packed_plus_minus(x, x_p, x_m, dim_1, dim_2, dim_full, offset)
!!
!!    Add double packed plus minus
!!    Written by Rolf H. Myhre, Eirik F. Kjønstad, and Sarai D. Folkestad Jan. 2021
!!
!!    x   :: packed 4d array
!!    x_p :: double packed plus combination 
!!    x_m :: double packed minus combination
!!
!!    x(pqrs) = x(pqrs) + x_p(qs,pr) + x_m(qs,pr)
!!    x(psrq) = x(psrq) + x_p(qs,pr) - x_m(qs,pr)
!!
      implicit none
!
      integer, intent(in) :: dim_1, dim_2, dim_full, offset
!
      complex(dp), dimension((dim_full*dim_1)*(dim_full*dim_1+1)/2), intent(inout) :: x
!
      complex(dp), dimension(dim_1*(dim_1+1)/2, dim_2*(dim_2+1)/2), intent(in) :: x_p
      complex(dp), dimension(dim_1*(dim_1-1)/2, dim_2*(dim_2-1)/2), intent(in) :: x_m
!
      integer :: p, r, q, s
      integer :: pr_p, pr_m, qs_p, qs_m
      integer :: pq, rs, ps, rq, pqrs, psrq
!
!$omp parallel do private(p,q,r,s,pq,rs,rq,ps,pqrs,psrq,pr_p,qs_p,pr_m,qs_m) &
!$omp schedule(guided)
      do s = 1, dim_1
         do q = 1, s
!
            qs_p = q + s*(s-1)/2
            qs_m = q + (s-1)*(s-2)/2
!
            do r = 1, dim_2
!
               rs = r + dim_full*(s-1) + offset
               rq = r + dim_full*(q-1) + offset
!
               do p = 1, r
!
                  pr_p = p + r*(r-1)/2
                  pr_m = p + (r-1)*(r-2)/2
!
                  pq = p + dim_full*(q-1) + offset
                  ps = p + dim_full*(s-1) + offset
!
                  pqrs = pq + rs*(rs-1)/2
                  psrq = max(rq,ps)*(max(rq,ps)-3)/2 + rq + ps
!
                  if (q .ne. s .and. p .ne. r) then
                     x(pqrs) = x(pqrs) + x_p(qs_p, pr_p) + x_m(qs_m, pr_m)
                     x(psrq) = x(psrq) + x_p(qs_p, pr_p) - x_m(qs_m, pr_m)
                  else
                     x(pqrs) = x(pqrs) + x_p(qs_p, pr_p)
                     x(psrq) = x(psrq) + x_p(qs_p, pr_p)
                  endif
!
               enddo
            enddo
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_double_packed_plus_minus
!
!
   subroutine add_single_packed_plus_minus(x, x_p, x_m, &
                                           dim_1, dim_2, dim_3, &
                                           dim_full, offset_2, offset_3)
!!
!!    Add single packed plus minus
!!    Written by Rolf H. Myhre, Eirik F. Kjønstad, and Sarai D. Folkestad Jan. 2021
!!
!!    x   :: packed 4d array
!!    x_p :: single packed plus combination
!!    x_m :: single packed minus combination
!!
!!    x(pqrs) = x(pqrs) + x_p(qs,p,r) + x_m(qs,p,r)
!!    x(psrq) = x(psrq) + x_p(qs,p,r) - x_m(qs,p,r)
!!
      implicit none
!
      integer, intent(in) :: dim_1, dim_2, dim_3, dim_full, offset_2, offset_3
!
      complex(dp), dimension((dim_full*dim_1)*(dim_full*dim_1+1)/2), intent(inout) :: x
!
      complex(dp), dimension(dim_1*(dim_1+1)/2, dim_2, dim_3), intent(in) :: x_p
      complex(dp), dimension(dim_1*(dim_1-1)/2, dim_2, dim_3), intent(in) :: x_m
!
      integer :: p, r, q, s, qs_p, qs_m
      integer :: pq, rs, pqrs, ps, rq, psrq
!
!$omp parallel do private(p,q,r,s,pq,rs,pqrs,ps,rq,psrq,qs_p,qs_m) &
!$omp schedule(guided)
      do s = 1, dim_1
         do q = 1, s
!
            qs_p = q + s*(s-1)/2
            qs_m = q + (s-1)*(s-2)/2
!
            do r = 1, dim_3
!
               rq = r + dim_full*(q-1) + offset_3
               rs = r + dim_full*(s-1) + offset_3
!
               do p = 1, dim_2
!
                  pq = p + dim_full*(q-1) + offset_2
                  ps = p + dim_full*(s-1) + offset_2
!
                  pqrs = max(pq,rs)*(max(pq,rs)-3)/2 + pq + rs
                  psrq = max(ps,rq)*(max(ps,rq)-3)/2 + ps + rq
!
                  if(q .ne. s) then
                     x(pqrs) = x(pqrs) + x_p(qs_p, p, r) + x_m(qs_m, p, r)
                     x(psrq) = x(psrq) + x_p(qs_p, p, r) - x_m(qs_m, p, r)
                  else
                     x(pqrs) = x(pqrs) + x_p(qs_p, p, r)
                     x(psrq) = x(psrq) + x_p(qs_p, p, r)
                  endif
!
               enddo
            enddo
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_single_packed_plus_minus
!
!
end module packed_array_utilities_c
