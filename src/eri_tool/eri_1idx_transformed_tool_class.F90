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
module eri_1idx_transformed_tool_class
!
!!
!! One-index transformed ERI tool class
!! Written by Sarai D. Folkestad, Sep 2021
!!
!! Handles the construction and delivery of the
!! one-index transformed electron repulsion integrals
!!
!! ERIs:
!!
!!    g_pqrs = L_pq^J L'_rs^J + L'_pq^J L_rs^J,
!!
!! where L' are transformed with respect to L
!!
!
   use parameters
!
   use global_out,                  only : output
   use memory_manager_class,        only : mem
   use abstract_eri_cholesky_class, only: abstract_eri_cholesky
   use range_class,                 only: range_
   use abstract_eri_tool_class,     only: abstract_eri_tool
!
   implicit none
!
   type, extends(abstract_eri_tool) :: eri_1idx_transformed_tool
!
      class(abstract_eri_cholesky), pointer :: L_1idx
!
   contains
!
      procedure, public :: get &
                        => get_eri_1idx_transformed_tool
!
      procedure, public :: update &
                        => update_eri_1idx_transformed_tool
!
      procedure, public :: get_memory_estimate &
                        => get_memory_estimate_eri_1idx_transformed_tool
!
      procedure, public :: get_memory_estimate_packed &
                        => get_memory_estimate_packed_eri_1idx_transformed_tool
!
      procedure, public :: get_symmetric_packed &
                        => get_symmetric_packed_eri_1idx_transformed_tool
!
   end type eri_1idx_transformed_tool
!
   interface eri_1idx_transformed_tool
!
      procedure :: new_eri_1idx_transformed_tool
!
   end interface eri_1idx_transformed_tool
!
contains
!
!
   function new_eri_1idx_transformed_tool(L, L_1idx) result(this)
!!
!!    New
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      type(eri_1idx_transformed_tool) :: this
      class(abstract_eri_cholesky), intent(in), target :: L, L_1idx
!
      this%L => L
      this%L_1idx => L_1idx
!
   end function new_eri_1idx_transformed_tool
!
!
   subroutine get_eri_1idx_transformed_tool(this, g_pqrs,    &
                                  first_p, last_p, &
                                  first_q, last_q, &
                                  first_r, last_r, &
                                  first_s, last_s, &
                                  alpha, beta)
!!
!!    Get
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(eri_1idx_transformed_tool), intent(inout) :: this
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
!
      real(dp), dimension(:,:,:), allocatable, target :: L_Jpq, L_Jrs
!
      call this%set_alpha_and_beta(alpha, beta, alpha_, beta_)
!
      call mem%alloc(L_Jpq, this%L%n_J, last_p-first_p+1, last_q-first_q+1)
      call this%L%get(L_Jpq, first_p, last_p, first_q, last_q)
!
      call mem%alloc(L_Jrs, this%L%n_J, last_r-first_r+1, last_s-first_s+1)
      call this%L_1idx%get(L_Jrs, first_r, last_r, first_s, last_s)
!
      call dgemm('T', 'N',                              &
                 (last_p-first_p+1)*(last_q-first_q+1), &
                 (last_r-first_r+1)*(last_s-first_s+1), &
                 this%L%n_J,                            &
                 alpha_,                                &
                 L_Jpq,                                 &
                 this%L%n_J,                            &
                 L_Jrs,                                 &
                 this%L%n_J,                            &
                 beta_,                                 &
                 g_pqrs,                                &
                 (last_p-first_p+1)*(last_q-first_q+1))
!
      call this%L_1idx%get(L_Jpq, first_p, last_p, first_q, last_q)
      call this%L%get(L_Jrs, first_r, last_r, first_s, last_s)
!
      call dgemm('T', 'N',                              &
                 (last_p-first_p+1)*(last_q-first_q+1), &
                 (last_r-first_r+1)*(last_s-first_s+1), &
                 this%L%n_J,                            &
                 alpha_,                                &
                 L_Jpq,                                 &
                 this%L%n_J,                            &
                 L_Jrs,                                 &
                 this%L%n_J,                            &
                 one,                                   &
                 g_pqrs,                                &
                 (last_p-first_p+1)*(last_q-first_q+1))
!
      call mem%dealloc(L_Jrs, this%L%n_J, last_r-first_r+1, last_s-first_s+1)
      call mem%dealloc(L_Jpq, this%L%n_J, last_p-first_p+1, last_q-first_q+1)
!
   end subroutine get_eri_1idx_transformed_tool
!
!
   subroutine get_symmetric_packed_eri_1idx_transformed_tool(this, g_pqpq,             &
                                                                    first_p, last_p,   &
                                                                    first_q, last_q,   &
                                                                    alpha, beta)
!!
!!    Get symmetric packed
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!
      use rectangular_full_packed_utilities_r, only: add_to_rfp
!
      implicit none
!
      class(eri_1idx_transformed_tool), intent(inout) :: this
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      real(dp), dimension((last_p - first_p + 1)            &
                           *(last_q - first_q + 1)          &
                           + mod((last_p - first_p + 1)     &
                           *(last_q - first_q + 1)+1, 2),   &
                          ((last_p - first_p + 1)           &
                           *(last_q - first_q + 1)+1)/2), intent(inout) :: g_pqpq
!
      real(dp), optional, intent(in) :: alpha, beta
!
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs
!
      call mem%alloc(g_pqrs, last_p-first_p+1, last_q-first_q+1, last_p-first_p+1, last_q-first_q+1)
!
      call this%get(g_pqrs,            &
                    first_p, last_p,   &
                    first_q, last_q,   &
                    first_p, last_p,   &
                    first_q, last_q,   &
                    alpha, beta)
!
      call add_to_rfp(one, g_pqrs, g_pqpq, (last_p-first_p+1)*(last_q-first_q+1))
      call mem%dealloc(g_pqrs, last_p-first_p+1, last_q-first_q+1, last_p-first_p+1, last_q-first_q+1)
!
   end subroutine get_symmetric_packed_eri_1idx_transformed_tool
!
!
   subroutine update_eri_1idx_transformed_tool(this)
!!
!!    Update
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(eri_1idx_transformed_tool), intent(inout) :: this
!
      call do_nothing(this)
!
   end subroutine update_eri_1idx_transformed_tool
!
!
   function get_memory_estimate_eri_1idx_transformed_tool(this,            &
                                                         first_p, last_p,  &
                                                         first_q, last_q,  &
                                                         first_r, last_r,  &
                                                         first_s, last_s) result(memory)
!!
!!    Get memory estimate
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(eri_1idx_transformed_tool), intent(in) :: this
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
      integer, intent(in) :: first_r, last_r
      integer, intent(in) :: first_s, last_s
!
      integer, dimension(2) :: memory
!
      integer :: dim_p, dim_q, dim_r, dim_s
!
      dim_p = last_p - first_p + 1
      dim_q = last_q - first_q + 1
      dim_r = last_r - first_r + 1
      dim_s = last_s - first_s + 1
!
      memory(1) = dim_p*dim_q*this%L%n_J + max(this%L%get_memory_estimate(first_p, last_p,       &
                                                                 first_q, last_q),               &
                                               this%L_1idx%get_memory_estimate(first_p, last_p,  &
                                                                           first_q, last_q))
      memory(2) = dim_r*dim_s*this%L%n_J + max(this%L%get_memory_estimate(first_r, last_r,       &
                                                                      first_s, last_s),          &
                                               this%L_1idx%get_memory_estimate(first_r, last_r,  &
                                                                      first_s, last_s))
!
   end function get_memory_estimate_eri_1idx_transformed_tool
!
!
   function get_memory_estimate_packed_eri_1idx_transformed_tool(this,  &
                                             first_p, last_p,           &
                                             first_q, last_q) result(memory)
!!
!!    Get memory estimate packed
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(eri_1idx_transformed_tool), intent(in) :: this
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      integer :: memory
!
      integer :: dim_p, dim_q
!
      dim_p = last_p - first_p + 1
      dim_q = last_q - first_q + 1
!
      memory = dim_p*dim_q*this%L%n_J + max(this%L%get_memory_estimate(first_p, last_p,         &
                                                                  first_q, last_q),             &
                                               this%L_1idx%get_memory_estimate(first_p, last_p, &
                                                                           first_q, last_q))
!
   end function get_memory_estimate_packed_eri_1idx_transformed_tool
!
!
end module eri_1idx_transformed_tool_class
