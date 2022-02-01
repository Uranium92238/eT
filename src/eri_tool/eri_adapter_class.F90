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
module eri_adapter_class
!
!!
!! ERI adapter class
!! Written by Sarai D. Folkestad, Aug-Sep 2021
!!
!! Used for simple integral block extraction
!! in CC code
!!
!! Acts as a wrapper for the eri tools
!!
!
   use parameters
!
   use global_out,                     only : output
   use memory_manager_class,           only : mem
   use abstract_eri_tool_class,        only: abstract_eri_tool
!
   implicit none
!
   type :: eri_adapter
!
      class(abstract_eri_tool), allocatable :: eri
!
      integer, private :: n_o, n_v, n_mo
      integer :: n_J
!
   contains
!
      procedure, public :: get &
                        => get_eri_adapter
!
      procedure, public :: get_packed &
                        => get_packed_eri_adapter
!
      procedure, public :: get_memory_estimate &
                        => get_memory_estimate_eri_adapter
!
      procedure, public :: get_memory_estimate_packed &
                        => get_memory_estimate_packed_eri_adapter
!
      procedure, private :: request_g
      procedure, private :: request_g_packed
      procedure, private :: index_setup
!
      final :: cleanup
!
   end type eri_adapter
!
   interface eri_adapter
!
      procedure :: new_eri_adapter
!
   end interface eri_adapter
!
contains
!
!
   function new_eri_adapter(eri, n_o, n_v) result(this)
!!
!!    New
!!    Written by Sarai D. Folkestad, Aug 2021
!!
      implicit none
!
      type(eri_adapter) :: this
      class(abstract_eri_tool), intent(in) :: eri
!
      integer :: n_o, n_v
!
      this%eri = eri
      this%n_o = n_o
      this%n_v = n_v
      this%n_mo = this%n_v + this%n_o
!
      this%n_J = this%eri%L%n_J
!
   end function new_eri_adapter
!
!
   pure subroutine index_setup(eri, x, first_index, last_index, &
                               first_restricted, last_restricted)
!!
!!    Index setup
!!    written by Rolf H. Myhre, Sep 2020
!!
!!    Sets first and last index based on the character x
!!
!!    x: string of characters to determine type of index
!!       o: occupied
!!       v: virtual
!!       f: full space
!!
!!    first_restricted, last_restricted : optional indices relative to range set by x
!!
!!    example: x='v' and first_restricted = 1 results in first_index = n_o + 1
!!
      implicit none
!
      class(eri_adapter),           intent(in)  :: eri
      character(len=1),             intent(in)  :: x
      integer,                      intent(out) :: first_index, last_index
      integer,            optional, intent(in)  :: first_restricted, last_restricted
!
      integer :: offset, default_first, default_last
!
      if (x .eq. 'o') then
!
         offset = 0
         default_first = 1
         default_last  = eri%n_o
!
      elseif (x .eq. 'f') then
!
         offset = 0
         default_first = 1
         default_last  = eri%n_mo
!
      elseif (x .eq. 'v') then
!
         offset = eri%n_o
         default_first = eri%n_o + 1
         default_last  = eri%n_mo
!
      else
!
         offset = -eri%n_mo - 1
         default_first = -1
         default_last  = -1
!
      endif
!
      if(present(first_restricted)) then
         first_index = first_restricted + offset
      else
         first_index = default_first
      endif
!
      if(present(last_restricted)) then
         last_index = last_restricted + offset
      else
         last_index = default_last
      endif
!
   end subroutine index_setup
!
!
   subroutine get_eri_adapter(this, x, g, &
                                     first_p, last_p, first_q, last_q, &
                                     first_r, last_r, first_s, last_s, &
                                     alpha, beta)
!!
!!    Get
!!    Written by Rolf H. Myhre, Sarai D. Folkestad and Eirik F. Kjønstad, 2020-2021
!!
      implicit none
!
      class(eri_adapter), intent(inout) :: this
!
      character(len=4), intent(in) :: x
!
      integer, optional, intent(in) :: first_p, last_p
      integer, optional, intent(in) :: first_q, last_q
      integer, optional, intent(in) :: first_r, last_r
      integer, optional, intent(in) :: first_s, last_s
!
      real(dp), intent(inout), dimension(1) :: g
!
      real(dp), optional, intent(in) :: alpha, beta
!
      integer :: full_first_p, full_last_p, full_first_q, full_last_q
      integer :: full_first_r, full_last_r, full_first_s, full_last_s
!
      call this%index_setup(x(1:1), full_first_p, full_last_p, first_p, last_p)
      call this%index_setup(x(2:2), full_first_q, full_last_q, first_q, last_q)
      call this%index_setup(x(3:3), full_first_r, full_last_r, first_r, last_r)
      call this%index_setup(x(4:4), full_first_s, full_last_s, first_s, last_s)
!
      call this%request_g(g, full_first_p, full_last_p, &
                         full_first_q, full_last_q, &
                         full_first_r, full_last_r, &
                         full_first_s, full_last_s, alpha, beta)
!
   end subroutine get_eri_adapter
!
!
   subroutine request_g(this, g, &
                    first_p, last_p, first_q, last_q, &
                    first_r, last_r, first_s, last_s, &
                    alpha, beta)
!!
!!    Request g
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!
      use reordering, only: sort_1234_to_1243, sort_1234_to_2134, sort_1234_to_2143
!
      implicit none
!
      class(eri_adapter), intent(inout) :: this
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
      integer, intent(in) :: first_r, last_r
      integer, intent(in) :: first_s, last_s
!
      real(dp), intent(inout), dimension((last_p - first_p + 1)*&
                                         (last_q - first_q + 1)*&
                                         (last_r - first_r + 1)*&
                                         (last_s - first_s + 1)), target :: g
!
      real(dp), optional, intent(in) :: alpha, beta
!
      call this%eri%get(g, first_p, last_p, &
                               first_q, last_q, &
                               first_r, last_r, &
                               first_s, last_s, alpha, beta)
!
   end subroutine request_g
!
!
   subroutine get_packed_eri_adapter(this, x, g,  &
                                           first_p, last_p,   &
                                           first_q, last_q,   &
                                           alpha, beta)
!!
!!    Get packed
!!    Written by Rolf H. Myhre, Sarai D. Folkestad and Eirik F. Kjønstad, 2020-2021
!!
      implicit none
!
      class(eri_adapter), intent(inout) :: this
!
      character(len=2), intent(in) :: x
!
      integer, optional, intent(in) :: first_p, last_p
      integer, optional, intent(in) :: first_q, last_q
!
      real(dp), intent(inout), dimension(1) :: g
!
      real(dp), optional, intent(in) :: alpha, beta
!
      integer :: full_first_p, full_last_p, full_first_q, full_last_q
!
      call this%index_setup(x(1:1), full_first_p, full_last_p, first_p, last_p)
      call this%index_setup(x(2:2), full_first_q, full_last_q, first_q, last_q)
!
      call this%request_g_packed(g, full_first_p, full_last_p, &
                         full_first_q, full_last_q, alpha, beta)
!
   end subroutine get_packed_eri_adapter
!
!
!
   subroutine request_g_packed(this, g,         &
                         first_p, last_p, &
                         first_q, last_q, &
                         alpha, beta)
!!
!!    Request g packed
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!
      use rectangular_full_packed_utilities_r, only: rfp_1234_to_2143
!
      implicit none
!
      class(eri_adapter), intent(inout) :: this
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      real(dp), dimension((last_p - first_p + 1)           &
                          *(last_q - first_q + 1)          &
                          + mod((last_p - first_p + 1)     &
                          *(last_q - first_q + 1)+1, 2),   &
                          (last_p - first_p + 1)           &
                          *((last_q - first_q + 1)+1)/2), intent(inout) :: g
!
      real(dp), optional, intent(in) :: alpha, beta
!
      call this%eri%get_symmetric_packed(g,                 &
                                         first_p, last_p,   &
                                         first_q, last_q,   &
                                         alpha, beta)
!
   end subroutine request_g_packed
!
!
   function get_memory_estimate_eri_adapter(this, x, &
                                          dim_p, dim_q, dim_r, dim_s) result(memory)
!!
!!    Get memory estimate
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Returns the memory per index pair pq and rs
!!    to get the integral g_pqrs
!!
!!    x: string of characters to determine type of index
!!       o: occupied
!!       v: virtual
!!       f: full space
!!
!!    p is a batching index: dim_p = 1
!!    p is not a batching index: dim_p given by the
!!                               corresponding character
!!                               in x
!!
!!
      implicit none
!
      class(eri_adapter), intent(in) :: this
!
      character(len=4), intent(in) :: x
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
      integer :: first_p, last_p
      integer :: first_q, last_q
      integer :: first_r, last_r
      integer :: first_s, last_s
!
      integer, dimension(2) :: memory
!
      call this%index_setup(x(1:1), first_p, last_p, 1, dim_p)
      call this%index_setup(x(2:2), first_q, last_q, 1, dim_q)
      call this%index_setup(x(3:3), first_r, last_r, 1, dim_r)
      call this%index_setup(x(4:4), first_s, last_s, 1, dim_s)
!
      memory = this%eri%get_memory_estimate(first_p, last_p, first_q, last_q, &
                                          first_r, last_r, first_s, last_s)
!
   end function get_memory_estimate_eri_adapter
!
!
   function get_memory_estimate_packed_eri_adapter(this, x, &
                                          dim_p, dim_q) result(memory)
!!
!!    Get memory estimate packed
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Returns the memory for index pair pq
!!    to get the packed symmetric integral g_pqrs
!!
!!    x: string of characters to determine type of index
!!       o: occupied
!!       v: virtual
!!       f: full space
!!
!!    p is a batching index: dim_p = 1
!!    p is not a batching index: dim_p given by the
!!                               corresponding character
!!                               in x
!!
!!
      implicit none
!
      class(eri_adapter), intent(in) :: this
!
      character(len=2), intent(in) :: x
!
      integer, intent(in) :: dim_p, dim_q
!
      integer :: first_p, last_p
      integer :: first_q, last_q
!
      integer :: memory
!
      call this%index_setup(x(1:1), first_p, last_p, 1, dim_p)
      call this%index_setup(x(2:2), first_q, last_q, 1, dim_q)
!
      memory = this%eri%get_memory_estimate_packed(first_p, last_p, first_q, last_q)
!
   end function get_memory_estimate_packed_eri_adapter
!
!
   subroutine cleanup(this)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, Jan 2022
!!
      implicit none
!
      type(eri_adapter), intent(inout) :: this
!
      if (allocated(this%eri)) deallocate(this%eri)
!
   end subroutine cleanup
!
!
end module eri_adapter_class
