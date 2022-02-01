!
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
module abstract_eri_tool_class
!
!!
!! Abstract ERI tool class
!! Writtan by Sarai D. Folkestad, Sep 2021
!!
!! Interface definition for the ERI tools
!!
!!
!
   use parameters
!
   use abstract_eri_cholesky_class, only: abstract_eri_cholesky
   use observer_class, only: observer
!
   implicit none
!
   type, extends(observer), abstract :: abstract_eri_tool
!
      class(abstract_eri_cholesky), pointer :: L
!
   contains
!
      procedure (get_abstract),                         deferred, public :: get
      procedure (get_symmetric_packed_abstract),        deferred, public :: get_symmetric_packed
      procedure (get_memory_estimate_abstract),         deferred, public :: get_memory_estimate
      procedure (get_memory_estimate_packed_abstract),  deferred, public :: get_memory_estimate_packed
!
      procedure, nopass :: set_alpha_and_beta
!
   end type abstract_eri_tool
!
   abstract interface
!
      subroutine get_abstract(this, g_pqrs,    &
                              first_p, last_p, &
                              first_q, last_q, &
                              first_r, last_r, &
                              first_s, last_s, &
                              alpha, beta)
!
         use parameters
         import abstract_eri_tool
!
         implicit none
!
         class(abstract_eri_tool), intent(inout) :: this
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
      end subroutine get_abstract
!
      subroutine get_symmetric_packed_abstract(this, g_pqpq,      &
                                               first_p, last_p,   &
                                               first_q, last_q,   &
                                               alpha, beta)
!
         use parameters
         import abstract_eri_tool
!
         implicit none
!
         class(abstract_eri_tool), intent(inout) :: this
!
         integer, intent(in) :: first_q, last_q
         integer, intent(in) :: first_p, last_p
!
!        g_pqpq in rectangular full packed (RFP) format.
         real(dp), dimension((last_p - first_p + 1)            &
                              *(last_q - first_q + 1)          &
                              + mod((last_p - first_p + 1)     &
                              *(last_q - first_q + 1)+1, 2),   &
                             ((last_p - first_p + 1)           &
                              *(last_q - first_q + 1)+1)/2), intent(inout) :: g_pqpq
!
         real(dp), optional, intent(in) :: alpha, beta

      end subroutine get_symmetric_packed_abstract
!
      subroutine update_abstract(this)
!
         import abstract_eri_tool
!
         implicit none
!
         class(abstract_eri_tool), intent(inout) :: this
!
      end subroutine update_abstract
!
      function get_memory_estimate_abstract(this,  &
                                                first_p, last_p, &
                                                first_q, last_q, &
                                                first_r, last_r, &
                                                first_s, last_s) result(memory)
!
         import abstract_eri_tool
!
         implicit none
!
         class(abstract_eri_tool), intent(in) :: this
!
         integer, intent(in) :: first_p, last_p
         integer, intent(in) :: first_q, last_q
         integer, intent(in) :: first_r, last_r
         integer, intent(in) :: first_s, last_s
!
         integer, dimension(2) :: memory
!
      end function get_memory_estimate_abstract
!
!
      function get_memory_estimate_packed_abstract(this,  &
                                                first_p, last_p, &
                                                first_q, last_q) result(memory)
!
         import abstract_eri_tool
!
         implicit none
!
         class(abstract_eri_tool), intent(in) :: this
!
         integer, intent(in) :: first_p, last_p
         integer, intent(in) :: first_q, last_q
!
         integer :: memory
!
      end function get_memory_estimate_packed_abstract
!
   end interface
!
contains
!
   subroutine set_alpha_and_beta(alpha, beta, alpha_, beta_)
!!
!!    Set alpha and beta
!!    Written by Sarai D. Folkestad, Oct 2021
!!
      implicit none
!
      real(dp), intent(in), optional :: alpha, beta
      real(dp), intent(out)          :: alpha_, beta_
!
      beta_ = zero
      if (present(beta)) beta_ = beta
!
      alpha_ = one
      if (present(alpha)) alpha_ = alpha
!
   end subroutine set_alpha_and_beta
!
end module abstract_eri_tool_class
