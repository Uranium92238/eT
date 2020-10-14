!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
!!    abstract eri tool class module
!!
!!    Written by Rolf H. Myhre, Sep 2020
!!
!!    Abstract tool with some general routines for real and complex eri tools
!!
!
   use parameters
   use global_out,           only : output
   use global_in,            only : input
   use memory_manager_class, only : mem
!
   implicit none
!
!
   type, abstract :: abstract_eri_tool
!
!     The number of Cholesky vectors, occupied, virtual and total mos
!
      integer :: n_J
      integer :: n_o
      integer :: n_v
      integer :: n_mo
!
!     Keep Cholesky vectors in memory
!
      logical :: cholesky_mem
      logical :: mo_eri_mem
!
   contains
!
      procedure :: read_general_settings => read_general_settings_abstract_eri_tool
!
!     Get ERI mem
!
      procedure :: get_eri_mem           => get_eri_mem_abstract_eri_tool
!
!     Various routines
!
      procedure :: room_for_cholesky     => room_for_cholesky_abstract_eri_tool
      procedure :: room_for_g_pqrs       => room_for_g_pqrs_abstract_eri_tool
!
      procedure :: is_pq_in_block        => is_pq_in_block_abstract_eri_tool
      procedure :: is_pq_contiguous      => is_pq_contiguous_abstract_eri_tool
!
      procedure :: index_setup           => index_setup_abstract_eri_tool
!
   end type abstract_eri_tool
!
!
contains
!
!
   subroutine read_general_settings_abstract_eri_tool(eri, eri_mem)
!!
!!    general settings
!!    Written by Eirik F. Kjønstad, Dec 2019
!!
      implicit none
!
      class(abstract_eri_tool) :: eri
!
      logical, intent(out) :: eri_mem
!
      character(len=200) :: cholesky_storage
      character(len=200) :: eri_storage
!
      if (input%requested_keyword_in_section('cholesky storage', 'integrals')) then
!
         call input%get_keyword_in_section('cholesky storage', &
                                           'integrals',        &
                                           cholesky_storage)
!
         if (cholesky_storage == 'memory') then
!
            eri%cholesky_mem = .true.
!
         elseif (cholesky_storage == 'disk') then
!
            eri%cholesky_mem = .false.
!
         else
!
            call output%error_msg('Did not recognize keyword value ' // trim(cholesky_storage) // &
                                  ' for cholesky storage.')
!
         endif
!
      endif
!
      if (input%requested_keyword_in_section('eri storage', 'integrals')) then
!
         call input%get_keyword_in_section('eri storage',     &
                                           'integrals',       &
                                            eri_storage)
!
        if (eri_storage == 'memory') then
!
            eri_mem = .true.
!
         elseif (eri_storage == 'none') then
!
            eri_mem = .false.
!
         else
!
            call output%error_msg('Did not recognize keyword value ' // trim(eri_storage) // &
                                  ' for eri storage.')
!
         endif
!
      endif
!
      if (input%requested_keyword_in_section('mo eri in memory', 'integrals')) then
!
            eri%mo_eri_mem = .true.
!
      endif
!
   end subroutine read_general_settings_abstract_eri_tool
!
!
   pure subroutine get_eri_mem_abstract_eri_tool(eri, string, req_pq, req_rs, &
                                                 dim_p, dim_q, dim_r, dim_s,  &
                                                 qp, sr)
!!
!!    Get ERI mem
!!    Written by Rolf H. Myhre Jun 2020
!!
!!    Adds the memory usage for get ERI depending on pq and rs to req_pq and req_rs
!!    Note that req_pq and req_rs must be initialized outside
!!    This routine is meant for estimates for the batching manager.
!!    If you are batching over index r, dim_r will typically be 1
!!    and the other dimensions should be full.
!!
!!    Temporary arrays needed for reordering, triggered by optional qp and sr
!!    are not allocated simultaneously in construct_eri, so this routine may
!!    overestimate total memory usage
!!
!!    See get_eri_mo for additional documentation
!!
      implicit none
!
      class(abstract_eri_tool), intent(in) :: eri
!
      character(len=4), intent(in) :: string
!
      integer, intent(inout) :: req_pq, req_rs
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      logical, optional, intent(in) :: qp, sr
!
      logical :: full_p, full_r
!
      full_p = ((string(1:1) .eq. 'v' .and. dim_p .eq. eri%n_v) .or. &
                (string(1:1) .eq. 'o' .and. dim_p .eq. eri%n_o))
!
      full_r = ((string(3:3) .eq. 'v' .and. dim_r .eq. eri%n_v) .or. &
                (string(3:3) .eq. 'o' .and. dim_r .eq. eri%n_o))
!
!     Can we use vectors in mem directly
      if(.not. (eri%cholesky_mem .and. full_p)) req_pq = req_pq + eri%n_J*dim_p*dim_q
      if(.not. (eri%cholesky_mem .and. full_r)) req_rs = req_rs + eri%n_J*dim_r*dim_s
!
!     Do we need temporary arrays for reordering
      if(present(qp)) then
         if(qp) req_pq = req_pq + eri%n_J*dim_p*dim_q
      endif
      if(present(sr)) then
         if(sr) req_rs = req_rs + eri%n_J*dim_r*dim_s
      endif
!
   end subroutine get_eri_mem_abstract_eri_tool
!
!
   pure function room_for_cholesky_abstract_eri_tool(eri, float_size) result(is_room)
!!
!!    Room for Cholesky
!!    Written by Eirik F. Kjønstad, Dec 2019
!!
!!    This routine is called to check whether the Cholesky vectors can be held in
!!    memory safely (< 20% of total available). If this is the case, the
!!    manager will keep a copy of them in memory.
!!
      implicit none
!
      class(abstract_eri_tool), intent(in) :: eri
!
      integer, intent(in) :: float_size
!
      integer(i64) :: required_mem
!
      integer, parameter :: fraction_of_total_mem = 5
!
      logical :: is_room
!
      is_room = .false.
!
      required_mem = int((eri%n_J)*(eri%n_mo)**2*float_size, kind=i64)
!
      if (required_mem .lt. mem%get_available()/fraction_of_total_mem) is_room = .true.
!
   end function room_for_cholesky_abstract_eri_tool
!
!
   pure function room_for_g_pqrs_abstract_eri_tool(eri, float_size) result(is_room)
!!
!!    Room for g_pqrs
!!    Written by Eirik F. Kjønstad, Jan 2019
!!
!!    This routine is called to check whether the T1-ERIs can be held in
!!    memory safely (< 20% of total available). If this is the case, the
!!    tool will keep a copy of g_pqrs in memory.
!!
      implicit none
!
      class(abstract_eri_tool), intent(in) :: eri
!
      integer, intent(in) :: float_size
!
      integer(i64) :: required_mem
!
      integer, parameter :: fraction_of_total_mem = 5
!
      logical :: is_room
!
      is_room = .false.
!
      required_mem = int((eri%n_mo**4)*float_size, kind=i64)
!
      if (required_mem .lt. mem%get_available()/fraction_of_total_mem) is_room = .true.
!
   end function room_for_g_pqrs_abstract_eri_tool
!
!
   pure function is_pq_in_block_abstract_eri_tool(eri, first_p, last_p, first_q, last_q) &
                 result(is_in_block)
!!
!!    is pq in block
!!    written by Rolf H. Myhre, Jun 2020
!!
!!    Check if all pq indices belong in a single block oo, vo, ov or vv
!!
      implicit none
!
      class(abstract_eri_tool), intent(in) :: eri
      integer, intent(in) :: first_p, last_p, first_q, last_q
!
      logical :: is_in_block
!
      if ((last_p .le. eri%n_o .or. first_p .gt. eri%n_o) .and. &
          (last_q .le. eri%n_o .or. first_q .gt. eri%n_o)) then
         is_in_block = .true.
      else
         is_in_block = .false.
      endif
!
   end function is_pq_in_block_abstract_eri_tool
!
!
   pure function is_pq_contiguous_abstract_eri_tool(eri, first_p, last_p, first_q, last_q) &
                 result(is_contiguous)
!!
!!    is pq contiguous
!!    written by Rolf H. Myhre, Jun 2020
!!
!!    Check if all pq in block and contiguous, i.e. full p index
!!
      implicit none
!
      class(abstract_eri_tool), intent(in) :: eri
      integer, intent(in) :: first_p, last_p, first_q, last_q
!
      logical :: is_contiguous
!
      if(eri%is_pq_in_block(first_p, last_p, first_q, last_q) .and. &
        (first_p .eq. 1 .or. first_p .eq. eri%n_o+1)          .and. &
        (last_p .eq. eri%n_o .or. last_p .eq. eri%n_mo)) then
         is_contiguous = .true.
      else
         is_contiguous = .false.
      endif
!
   end function is_pq_contiguous_abstract_eri_tool
!
!
   pure subroutine index_setup_abstract_eri_tool(eri, x, first_index, last_index, opt_first, opt_last)
!!
!!    index setup
!!    written by Rolf H. Myhre, Sep 2020
!!
!!    Sets first and last index based on the character x
!!
!!    x: character to determine type of index
!!       o: occupied
!!       v: virtual
!!       f: full space
!!
!!    opt_first, opt_last : optional indices relative to range set by x
!!
!!    example: x='v' and opt_first = 1 results in first_index = n_o + 1
!!
      implicit none
!
      class(abstract_eri_tool), intent(in) :: eri
!
      character(len=1), intent(in) :: x
!
      integer, intent(out) :: first_index, last_index
      integer, optional, intent(in) :: opt_first, opt_last
!
      integer :: offset, default_first, default_last
!
      if (x .eq. 'o') then !occupied
         offset = 0
         default_first = 1
         default_last  = eri%n_o
      elseif (x .eq. 'f') then !full
         offset = 0
         default_first = 1
         default_last  = eri%n_mo
      elseif (x .eq. 'v') then !virtual
         offset = eri%n_o
         default_first = eri%n_o + 1
         default_last  = eri%n_mo
      else
         offset = -eri%n_mo - 1
         default_first = -1
         default_last  = -1
      endif
!
      if(present(opt_first)) then
         first_index = opt_first + offset
      else
         first_index = default_first
      endif
!
      if(present(opt_last)) then
         last_index = opt_last + offset
      else
         last_index = default_last
      endif
!
   end subroutine index_setup_abstract_eri_tool
!
!
end module abstract_eri_tool_class
