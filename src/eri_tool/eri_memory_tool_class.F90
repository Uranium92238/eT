!
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
module eri_memory_tool_class
!
!!
!! ERI memory tool
!! Written by Sarai D. Folkestad, Sep 2021
!!
!! Handles the memory storage and delivery of the
!! electron repulsion integrals
!!
!! ERIs:
!!
!!    g_pqrs = L_pq^J L_rs^J
!!
!!
!
   use parameters
!
   use global_out,                   only : output
   use memory_manager_class,         only : mem
   use abstract_eri_cholesky_class,  only: abstract_eri_cholesky
   use range_class,                  only: range_
   use abstract_eri_tool_class,      only: abstract_eri_tool
!
   implicit none
!
   type, extends(abstract_eri_tool) :: eri_memory_tool
!
      real(dp), dimension(:,:,:,:), allocatable :: g
!
   contains
!
      procedure, public :: get &
                        => get_eri_memory_tool
!
      procedure, public :: get_symmetric_packed &
                        => get_symmetric_packed_eri_memory_tool
!
      procedure, public :: update &
                        => update_eri_memory_tool
!
      procedure, public :: get_memory_estimate &
                        => get_memory_estimate_eri_memory_tool
!
      procedure, public :: get_memory_estimate_packed &
                        => get_memory_estimate_packed_eri_memory_tool
!
      final :: cleanup
!
   end type eri_memory_tool
!
   interface eri_memory_tool
!
      procedure :: new_eri_memory_tool
!
   end interface eri_memory_tool
!
contains
!
   function new_eri_memory_tool(L) result(this)
!!
!!    New
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      type(eri_memory_tool) :: this
      class(abstract_eri_cholesky), intent(in), target :: L
!
      this%L => L
!
   end function new_eri_memory_tool
!
!
   subroutine get_eri_memory_tool(this, g_pqrs,    &
                                         first_p, last_p, &
                                         first_q, last_q, &
                                         first_r, last_r, &
                                         first_s, last_s, &
                                         alpha, beta)
!
!!
!!    Get
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!       g_pqrs = beta * g_pqrs + alpha * g
!!
!!    alpha is default one
!!    beta is default zero
!!
      implicit none
!
      class(eri_memory_tool), intent(inout) :: this
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
      integer, intent(in) :: first_r, last_r
      integer, intent(in) :: first_s, last_s
!
      real(dp), dimension(first_p:last_p, &
                          first_q:last_q, &
                          first_r:last_r, &
                          first_s:last_s), intent(inout) :: g_pqrs
!
      real(dp), optional, intent(in) :: alpha, beta
!
      real(dp) :: alpha_, beta_
      integer  :: p, q, r, s
!
      call this%set_alpha_and_beta(alpha, beta, alpha_, beta_)
!
      if (beta_ .eq. zero) then
!$omp parallel do private (s, r, q, p)
         do s = first_s, last_s
            do r = first_r, last_r
               do q = first_q, last_q
                  do p = first_p, last_p
                     g_pqrs(p, q, r, s) = alpha_*this%g(p, q, r, s )
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
      else
!$omp parallel do private (s, r, q, p)
         do s = first_s, last_s
            do r = first_r, last_r
               do q = first_q, last_q
                  do p = first_p, last_p
                     g_pqrs(p, q, r, s) = beta_*g_pqrs(p, q, r, s) + alpha_*this%g(p, q, r, s )
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
      endif
!
   end subroutine get_eri_memory_tool
!
!
   subroutine get_symmetric_packed_eri_memory_tool(this, g_pqpq,      &
                                                          first_p, last_p,   &
                                                          first_q, last_q,   &
                                                          alpha, beta)
!!
!!    Get symmetric packed
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    Returns the array in upper rectangular full packed format
!!
!!       g_pqpq = beta * g_pqpq + alpha * g
!!
!!    alpha is default one
!!    beta is default zero
!!
!!    Returns the symmetric integral in the rectangular full packed format
!!
      implicit none
!
      class(eri_memory_tool), intent(inout) :: this
!
      integer, intent(in) :: first_q, last_q
      integer, intent(in) :: first_p, last_p
!
      real(dp), dimension((last_p - first_p + 1)            &
                           *(last_q - first_q + 1)          &
                           + mod((last_p - first_p + 1)     &
                           *(last_q - first_q + 1)+1, 2),   &
                          ((last_p - first_p + 1)           &
                           *(last_q - first_q + 1)+1)/2), intent(inout) :: g_pqpq
!
      real(dp), optional, intent(in) :: alpha, beta
      real(dp) :: alpha_, beta_
!
      integer :: dim_p, dim_q, dim
      integer :: y, x, pq, rs, p, q, r, s
!
      dim_p = last_p - first_p + 1
      dim_q = last_q - first_q + 1
      dim = dim_p*dim_q
!
      call this%set_alpha_and_beta(alpha, beta, alpha_, beta_)
!
      if(beta_ .eq. zero) then
!
!$omp parallel do private(x, y, pq, rs, p, q, r, s)
         do y = 1, (dim_p*dim_q+1)/2
            do x = 1, dim_p*dim_q + mod(dim_p*dim_q+1,2)
!
               if (x .le. y + dim/2) then
                  pq = x
                  rs = y + dim/2
               else
                  pq = y
                  rs = x - dim/2 - 1
               endif
!
               p = mod(pq-1,dim_p)
               r = mod(rs-1,dim_p)
               q = (pq-1)/dim_p
               s = (rs-1)/dim_p
!
               g_pqpq(x,y) = alpha_*this%g(first_p+p, first_q+q, first_p+r, first_q+s)
!
            enddo
         enddo
!$omp end parallel do
!
      else
!
!$omp parallel do private(x, y, pq, rs, p, q, r, s)
         do y = 1, (dim_p*dim_q+1)/2
            do x = 1, dim_p*dim_q + mod(dim_p*dim_q+1,2)
!
               if (x .le. y + dim/2) then
                  pq = x
                  rs = y + dim/2
               else
                  pq = y
                  rs = x - dim/2 - 1
               endif
!
               p = mod(pq-1,dim_p)
               r = mod(rs-1,dim_p)
               q = (pq-1)/dim_p
               s = (rs-1)/dim_p
!
               g_pqpq(x,y) =  beta_*g_pqpq(x,y) &
                           + alpha_*this%g(first_p+p,first_q+q,first_p+r,first_q+s)
!
            enddo
         enddo
!$omp end parallel do
!
      endif
!
   end subroutine get_symmetric_packed_eri_memory_tool
!
!
   subroutine update_eri_memory_tool(this)
!!
!!    Update
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(eri_memory_tool), intent(inout) :: this
!
      real(dp), dimension(:,:,:), pointer :: L
!
      if (.not. allocated(this%g)) call mem%alloc(this%g, this%L%dim_, this%L%dim_, this%L%dim_, this%L%dim_)
!
      call this%L%load_block(L, 1, this%L%dim_, 1, this%L%dim_)
!
      call dgemm('T', 'N',       &
                 this%L%dim_**2, &
                 this%L%dim_**2, &
                 this%L%n_J,     &
                 one,            &
                 L,              &
                 this%L%n_J,     &
                 L,              &
                 this%L%n_J,     &
                 zero,           &
                 this%g,         &
                 this%L%dim_**2)
!
      call this%L%offload_block(1, this%L%dim_, 1, this%L%dim_)
!
   end subroutine update_eri_memory_tool
!
!
   subroutine cleanup(this)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      type(eri_memory_tool), intent(inout) :: this
!
      if (allocated(this%g)) call mem%dealloc(this%g, this%L%dim_, this%L%dim_, this%L%dim_, this%L%dim_)
!
   end subroutine cleanup
!
!
   function get_memory_estimate_eri_memory_tool(this,  &
                                             first_p, last_p, &
                                             first_q, last_q, &
                                             first_r, last_r, &
                                             first_s, last_s) result(memory)
!!
!!    Get memory estimate
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Given the first_x and last_q x, x = p,q,r,s
!!    asks the Cholesky tool if extra memory is needed
!!
!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(eri_memory_tool), intent(in) :: this
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
      integer, intent(in) :: first_r, last_r
      integer, intent(in) :: first_s, last_s
!
      integer, dimension(2) :: memory
!
      memory = 0
!
      call do_nothing(this)
      call do_nothing([first_p, last_p, &
                       first_q, last_q, &
                       first_r, last_r, &
                       first_s, last_s])
!
   end function get_memory_estimate_eri_memory_tool
!
!
   function get_memory_estimate_packed_eri_memory_tool(this,  &
                                             first_p, last_p, &
                                             first_q, last_q) result(memory)
!!
!!    Get memory estimate packed
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Given the first_x and last_q x, x = p,q
!!    asks the Cholesky tool if extra memory is needed
!!
!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(eri_memory_tool), intent(in) :: this
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      integer :: memory
!
      memory = 0
!
      call do_nothing(this)
      call do_nothing([first_p, last_p, &
                       first_q, last_q])
!
   end function get_memory_estimate_packed_eri_memory_tool
!
!
end module eri_memory_tool_class
