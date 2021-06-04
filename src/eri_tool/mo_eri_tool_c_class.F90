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
module mo_eri_tool_c_class
!
!!
!!    MO eri tool complex class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Extensively modified by Rolf H. Myhre, June 2020
!!
!!    Tool that handles the two-electron repulsion integrals in the MO/T1/"C1"
!!    basis. The tool uses Cholesky vectors L_pq^J to construct
!!    the electron repulsion integrals:
!!
!!       g_pqrs = sum_J L_pq^J L_rs^J
!!
!!    Calling get_eri_mo returns the integrals
!!    Calling get_eri_mo_mem returns an estimate for the amount of memory needed
!!
!
   use global_out,   only : output
   use global_in,    only : input
!
   use abstract_eri_tool_class,  only : abstract_eri_tool
   use mo_eri_tool_class,        only : mo_eri_tool
!
   use memory_manager_class,     only : mem
   use direct_stream_file_class, only : direct_stream_file
   use timings_class,            only : timings
   use batching_index_class,     only : batching_index
!
   use reordering
!
   implicit none
!
!
   type, extends(abstract_eri_tool) :: mo_eri_tool_c
!
!     Files to keep L_pq^J in MO and T1 bases
!
      type(direct_stream_file) :: cholesky_mo
!
!     Arrays for MO Cholesky vectors
!
      complex(dp), dimension(:,:,:), allocatable :: L_J_oo_mo
      complex(dp), dimension(:,:,:), allocatable :: L_J_vo_mo
      complex(dp), dimension(:,:,:), allocatable :: L_J_ov_mo
      complex(dp), dimension(:,:,:), allocatable :: L_J_vv_mo
!
!     Keep full electron repulsion matrix g_pqrs in memory? (+ allocatable arrays for keeping it)
!
      complex(dp), dimension(:,:,:,:), allocatable :: g_pqrs_mo
!
   contains
!
!     Read, print settings
!
      procedure :: print_settings            => print_settings_mo_eri_tool_c
      procedure :: read_settings             => read_settings_mo_eri_tool_c
!
!     Cleanup
!
      procedure :: cleanup                   => cleanup_mo_eri_tool_c
!
!     Initialization
!
      procedure :: initialize                => initialize_mo_eri_tool_c
      procedure :: copy_from_mo              => copy_from_mo_mo_eri_tool_c
      procedure :: copy_from_real_mo         => copy_from_real_mo_mo_eri_tool_c
!
!     write, read, copy_in, copy_out vectors
!
      procedure :: write_cholesky            => write_cholesky_mo_eri_tool_c
      procedure :: read_cholesky             => read_cholesky_mo_eri_tool_c
!
      procedure :: copy_in_cholesky          => copy_in_cholesky_mo_eri_tool_c
      procedure :: copy_out_cholesky         => copy_out_cholesky_mo_eri_tool_c
!
!     MO Cholesky veectors
!
      procedure :: set_cholesky_mo           => set_cholesky_mo_mo_eri_tool_c
      procedure :: get_cholesky_mo           => get_cholesky_mo_mo_eri_tool_c
!
      procedure :: update_cholesky_mo        => update_cholesky_mo_mo_eri_tool_c
!
!     Get and construct routines for g_pqrs in MO basis
!
      procedure :: get_eri_mo_mem            => get_eri_mo_mem_mo_eri_tool_c
      procedure :: get_eri_mo                => get_eri_mo_mo_eri_tool_c
      procedure :: get_g_pqrs_mo             => get_g_pqrs_mo_mo_eri_tool_c
      procedure :: construct_g_pqrs_mo       => construct_g_pqrs_mo_mo_eri_tool_c
!
      procedure :: copy_g_pqrs               => copy_g_pqrs_mo_eri_tool_c
      procedure :: copy_g_pqrs_to_packed     => copy_g_pqrs_to_packed_mo_eri_tool_c
!
!     Various routines
!
      procedure :: place_g_mo_in_memory      => place_g_mo_in_memory_mo_eri_tool_c
!
      procedure :: set_pointer_mo            => set_pointer_mo_mo_eri_tool_c
      procedure :: L_pointer_setup_mo        => L_pointer_setup_mo_mo_eri_tool_c
      procedure :: construct_g_from_L        => construct_g_from_L_mo_eri_tool_c
      procedure :: construct_g_symm_from_L   => construct_g_symm_from_L_mo_eri_tool_c
      procedure :: construct_g_packed_from_L => construct_g_packed_from_L_mo_eri_tool_c
!
   end type mo_eri_tool_c
!
!
   interface mo_eri_tool_c
!
      procedure :: new_mo_eri_tool_c
      procedure :: new_mo_eri_tool_c_from_real
      procedure :: new_mo_eri_tool_c_from_complex
!
   end interface mo_eri_tool_c
!
!
contains
!
!
   function new_mo_eri_tool_c(n_o, n_v, n_J, need_g_abcd) result(eri)
!!
!!    New MO integral tool
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Creates the MO integral tool and sets internal variables.
!!
!!    n_o:          Number of occupied orbitals
!!    n_v:          Number of virtual orbitals
!!    n_J:          Number of Cholesky vectors
!!    need_g_abcd:  Will we have use of the full v^4 integrals? If so,
!!                  default is to store the full ERI matrix in memory
!!                  if there is plenty of room.
!!
      implicit none
!
      type(mo_eri_tool_c) :: eri
!
      integer, intent(in) :: n_o
      integer, intent(in) :: n_v
      integer, intent(in) :: n_J
      logical, intent(in) :: need_g_abcd
!
      eri%n_J  = n_J
      eri%n_o  = n_o
      eri%n_v  = n_v
      eri%n_mo = n_o + n_v
!
      eri%cholesky_mo = direct_stream_file('cholesky_mo', eri%n_J)
!
!     Default: Cholesky in memory, integral matrix in memory
!              depends on amount of memory and "need_g_abcd"
!
      if (eri%room_for_cholesky(2*dp)) then
!
         eri%cholesky_mem = .true.
!
      else
!
         eri%cholesky_mem = .false.
!
      endif
!
      if (need_g_abcd .and. eri%room_for_g_pqrs(2*dp)) then
!
        eri%mo_eri_mem = .true.
!
      else
!
        eri%mo_eri_mem = .false.
!
      endif
!
!     Read non-default settings
!
      call eri%read_settings()
!
!     Print settings
!
      call eri%print_settings()
!
   end function new_mo_eri_tool_c
!
!
   function new_mo_eri_tool_c_from_complex(template) result(eri)
!!
!!    New MO integral tool complex from complex
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Creates the MO integral tool and sets internal variables based on template MO tool.
!!
      implicit none
!
      type(mo_eri_tool_c) :: eri
!
      class(mo_eri_tool_c) :: template
!
      eri%n_J  = template%n_J
      eri%n_o  = template%n_o
      eri%n_v  = template%n_v
      eri%n_mo = template%n_mo
!
      eri%cholesky_mo = direct_stream_file('cholesky_mo', eri%n_J)
!
!     Memory logicals - Cholesky in memory (or file)? ERI in memory (or file)?
!
      eri%cholesky_mem = template%cholesky_mem
      eri%mo_eri_mem      = template%mo_eri_mem
!
!     Print settings
!
      call eri%print_settings()
!
   end function new_mo_eri_tool_c_from_complex
!
!
   function new_mo_eri_tool_c_from_real(template) result(eri)
!!
!!    New MO integral tool complex from real
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Creates the integral tool based on template.
!!
      implicit none
!
      type(mo_eri_tool_c) :: eri
!
      class(mo_eri_tool) :: template
!
      eri%n_J  = template%n_J
      eri%n_o  = template%n_o
      eri%n_v  = template%n_v
      eri%n_mo = template%n_mo
!
      eri%cholesky_mo = direct_stream_file('cholesky_mo', eri%n_J)
!
!     Memory logicals - Cholesky in memory (or file)? ERI in memory (or file)?
!
      eri%cholesky_mem = template%cholesky_mem
      eri%mo_eri_mem      = template%mo_eri_mem
!
!     Print settings
!
      call eri%print_settings()
!
   end function new_mo_eri_tool_c_from_real
!
!
   subroutine print_settings_mo_eri_tool_c(eri)
!!
!!    Print settings
!!    Written by Eirik F. Kjønstad, Dec 2019
!!
      implicit none
!
      class(mo_eri_tool_c) :: eri
!
      call output%printf('m', '- Settings for integral handling:', fs='(/t3,a)')
!
      call output%printf('m', &
                         'Cholesky vectors in memory: (l0)', &
                         logs=[eri%cholesky_mem],            &
                         fs='(/t6,a)')
!
      call output%printf('m', &
                         'ERI matrix in memory:       (l0)', &
                         logs=[eri%mo_eri_mem],                 &
                         fs='(t6,a)')
!
   end subroutine print_settings_mo_eri_tool_c
!
!
   subroutine read_settings_mo_eri_tool_c(eri)
!!
!!    Read settings
!!    Written by Eirik F. Kjønstad, Dec 2019
!!
      implicit none
!
      class(mo_eri_tool_c) :: eri
!
      call eri%read_general_settings(eri%mo_eri_mem)
!
   end subroutine read_settings_mo_eri_tool_c
!
!
   subroutine initialize_mo_eri_tool_c(eri)
!!
!!    Initialize
!!    Written by Rolf H. Myhre, Aug 2020
!!
!!    Allocates arrays for Cholesky vectors and the ERIs
!!    if they are to be stored in memory.
!!
      implicit none
!
      class(mo_eri_tool_c) :: eri
!
!     Allocate Cholesky arrays
!
      if (eri%cholesky_mem) then
!
         call mem%alloc(eri%L_J_oo_mo, eri%n_J, eri%n_o, eri%n_o)
         call mem%alloc(eri%L_J_ov_mo, eri%n_J, eri%n_o, eri%n_v)
         call mem%alloc(eri%L_J_vo_mo, eri%n_J, eri%n_v, eri%n_o)
         call mem%alloc(eri%L_J_vv_mo, eri%n_J, eri%n_v, eri%n_v)
!
      endif
!
!        Allocate integrals
!
      if (eri%mo_eri_mem) then
         call mem%alloc(eri%g_pqrs_mo, eri%n_mo, eri%n_mo, &
                                       eri%n_mo, eri%n_mo)
      endif
!
   end subroutine initialize_mo_eri_tool_c
!
!
   subroutine copy_from_mo_mo_eri_tool_c(eri, template)
!!
!!    Copy from MO
!!    Written by Rolf H. Myhre, Aug 2020
!!
!!    Copies MO Cholesky and MO ERIs if in memory
!!
      implicit none
!
      class(mo_eri_tool_c) :: eri
!
      class(mo_eri_tool_c) :: template
!
!     Copy Cholesky vectors from template - if present
!
      if (eri%cholesky_mem .and. template%cholesky_mem) then
!
         call zcopy(eri%n_J*eri%n_o*eri%n_o, &
                    template%L_J_oo_mo,  &
                    1, eri%L_J_oo_mo, 1)
!
         call zcopy(eri%n_J*eri%n_o*eri%n_v, &
                    template%L_J_ov_mo,  &
                    1, eri%L_J_ov_mo, 1)
!
         call zcopy(eri%n_J*eri%n_v*eri%n_o, &
                    template%L_J_vo_mo,  &
                    1, eri%L_J_vo_mo, 1)
!
         call zcopy(eri%n_J*eri%n_v*eri%n_v, &
                    template%L_J_vv_mo,  &
                    1, eri%L_J_vv_mo, 1)
!
      endif
!
!     Copy integrals from template - if present
!
      if (eri%mo_eri_mem .and. template%mo_eri_mem) then
!
         call zcopy(eri%n_mo**4,        &
                    template%g_pqrs_mo, &
                    1,                  &
                    eri%g_pqrs_mo,      &
                    1)
!
      endif
!
   end subroutine copy_from_mo_mo_eri_tool_c
!
!
   subroutine copy_from_real_mo_mo_eri_tool_c(eri, template)
!!
!!    Copy from real MO
!!    Written by Rolf H. Myhre, Aug 2020
!!
!!    Copies MO Cholesky and MO ERIs if in memory
!!
      implicit none
!
      class(mo_eri_tool_c) :: eri
!
      class(mo_eri_tool)   :: template
!
!     Copy Cholesky vectors from template - if present
!
      if (eri%cholesky_mem .and. template%cholesky_mem) then
!
         eri%L_J_oo_mo = cmplx(template%L_J_oo_mo, zero, dp)
         eri%L_J_vo_mo = cmplx(template%L_J_vo_mo, zero, dp)
         eri%L_J_ov_mo = cmplx(template%L_J_ov_mo, zero, dp)
         eri%L_J_vv_mo = cmplx(template%L_J_vv_mo, zero, dp)
!
      endif
!
!     Copy integrals from template - if present
!
      if (eri%mo_eri_mem .and. template%mo_eri_mem) then
!
         eri%g_pqrs_mo = cmplx(template%g_pqrs_mo, zero, dp)
!
      endif
!
   end subroutine copy_from_real_mo_mo_eri_tool_c
!
!
   subroutine write_cholesky_mo_eri_tool_c(eri, L_J_pq, &
                                      first_p, last_p, first_q, last_q, &
                                      file_)
!!
!!    Write cholesky vectors
!!    Written by Rolf H. Myhre, May 2020
!!
!!    Writes general Cholesky vectors L_J_pq to file distributed into oo, vo, ov, and vv blocks
!!
!!    first_p, last_p: first and last full index of p
!!    first_q, last_q: first and last full index of q
!!
!!    file_: file to write to
!!
      implicit none
!
      class(mo_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      complex(dp), dimension(eri%n_J, first_p:last_p, first_q:last_q), intent(in) :: L_J_pq
!
      class(direct_stream_file) :: file_
!
      integer :: b, j, offset
!
      if (first_q .le. eri%n_o) then !write to oo and/or vo
!
         if (first_p .le. eri%n_o) then !write to oo
            offset = 0
            if(first_p .eq. 1 .and. last_p .eq. eri%n_o) then
!
               call file_%write_(L_J_pq(:, :, first_q:min(last_q, eri%n_o)), &
                                 offset + (first_q-1)*eri%n_o + 1,           &
                                 offset + last_q*eri%n_o)
            else
               do j = first_q, min(last_q, eri%n_o)
!
                  call file_%write_(L_J_pq(:, first_p:min(last_p, eri%n_o), j),      &
                                    offset + (j - 1)*eri%n_o + first_p,              &
                                    offset + (j - 1)*eri%n_o + min(last_p, eri%n_o))
               enddo
            endif
         endif
!
         if (last_p .gt. eri%n_o) then !write to vo
            offset = eri%n_o**2
            if(first_p .eq. eri%n_o + 1 .and. last_p .eq. eri%n_mo) then
!
               call file_%write_(L_J_pq(:, :, first_q:min(last_q, eri%n_o)), &
                                 offset + (first_q - 1)*eri%n_v + 1,         &
                                 offset + last_q*eri%n_v)
            else
               do j = first_q, min(last_q, eri%n_o)
!
                  call file_%write_(L_J_pq(:, max(first_p, eri%n_o + 1):last_p, j),       &
                                    offset + (j - 1)*eri%n_v + max(first_p - eri%n_o, 1), &
                                    offset + (j - 1)*eri%n_v + last_p - eri%n_o)
               enddo
            endif
         endif
      endif
!
      if (last_q .gt. eri%n_o) then !write to ov and/or vv
!
         if (first_p .le. eri%n_o) then !Write to ov
            offset = eri%n_o**2 + eri%n_o*eri%n_v
            if(first_p .eq. 1 .and. last_p .eq. eri%n_o) then
!
               call file_%write_(L_J_pq(:, :, max(first_q, eri%n_o+1):last_q), &
                                 offset + (first_q - eri%n_o - 1)*eri%n_o + 1, &
                                 offset + (last_q - eri%n_o)*eri%n_o)
            else
               do b = max(first_q, eri%n_o + 1), last_q
!
                  call file_%write_(L_J_pq(:, first_p:min(last_p, eri%n_o), b),                &
                                    offset + (b - eri%n_o - 1)*eri%n_o + first_p,              &
                                    offset + (b - eri%n_o - 1)*eri%n_o + min(last_p, eri%n_o))
               enddo
            endif
         endif
!
         if (last_p .gt. eri%n_o) then !Write to vv
            offset = eri%n_o**2 + 2*eri%n_o*eri%n_v
            if(first_p .eq. eri%n_o+1 .and. last_p .eq. eri%n_mo) then
!
               call file_%write_(L_J_pq(:, :, max(first_q, eri%n_o + 1):last_q), &
                                 offset + (first_q - eri%n_o - 1)*eri%n_v + 1,   &
                                 offset + (last_q - eri%n_o)*eri%n_v)
            else
               do b = max(first_q, eri%n_o + 1), last_q
!
                  call file_%write_(L_J_pq(:, max(first_p, eri%n_o + 1):last_p, b), &
                                    offset + (b-eri%n_o - 1)*eri%n_v + max(first_p - eri%n_o,1), &
                                    offset + (b-eri%n_o - 1)*eri%n_v + last_p - eri%n_o)
               enddo
            endif
         endif
      endif
!
   end subroutine write_cholesky_mo_eri_tool_c
!
!
   subroutine read_cholesky_mo_eri_tool_c(eri, L_J_pq, &
                                     first_p, last_p, first_q, last_q, &
                                     file_)
!!
!!    read cholesky vectors
!!    Written by Rolf H. Myhre, May 2020
!!
!!    reads general Cholesky vectors L_J_pq to file distributed into oo, vo, ov, and vv blocks
!!
!!    first_p, last_p: first and last full index of p
!!    first_q, last_q: first and last full index of q
!!
!!    file_: file to read from
!!
      implicit none
!
      class(mo_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      complex(dp), dimension(eri%n_J, first_p:last_p, first_q:last_q), intent(out) :: L_J_pq
!
      class(direct_stream_file) :: file_
!
      integer :: b, j, offset
!
      if (first_q .le. eri%n_o) then !Read from oo and/or vo
!
         if (first_p .le. eri%n_o) then !Read from oo
            offset = 0
!
            if(first_p .eq. 1 .and. last_p .eq. eri%n_o) then
!
               call file_%read_(L_J_pq(:, :, first_q:min(last_q, eri%n_o)), &
                                offset + (first_q - 1)*eri%n_o + 1,         &
                                offset + last_q*eri%n_o)
            else
               do j = first_q, min(last_q, eri%n_o)
!
                  call file_%read_(L_J_pq(:, first_p:min(last_p, eri%n_o), j),    &
                                   offset + (j-1)*eri%n_o + first_p,              &
                                   offset + (j-1)*eri%n_o + min(last_p, eri%n_o))
               enddo
            endif
         endif
!
         if (last_p .gt. eri%n_o) then !Read from vo
            offset = eri%n_o**2
            if(first_p .eq. eri%n_o+1 .and. last_p .eq. eri%n_mo) then
!
               call file_%read_(L_J_pq(:, :, first_q:min(last_q, eri%n_o)), &
                                offset + (first_q-1)*eri%n_v + 1,           &
                                offset + last_q*eri%n_v)
            else
               do j = first_q, min(last_q, eri%n_o)
!
                  call file_%read_(L_J_pq(:, max(first_p, eri%n_o + 1):last_p, j),     &
                                   offset + (j-1)*eri%n_v + max(first_p - eri%n_o, 1), &
                                   offset + (j-1)*eri%n_v + last_p - eri%n_o)
               enddo
            endif
         endif
      endif
!
      if (last_q .gt. eri%n_o) then !read from ov and/or vv
!
         if (first_p .le. eri%n_o) then !Read from ov
            offset = eri%n_o**2 + eri%n_o*eri%n_v
            if(first_p .eq. 1 .and. last_p .eq. eri%n_o) then
!
               call file_%read_(L_J_pq(:, :, max(first_q, eri%n_o+1):last_q), &
                                offset + (first_q - eri%n_o - 1)*eri%n_o + 1, &
                                offset + (last_q - eri%n_o)*eri%n_o)
            else
               do b = max(first_q, eri%n_o + 1), last_q
!
                  call file_%read_(L_J_pq(:, first_p:min(last_p, eri%n_o), b),              &
                                   offset + (b-eri%n_o - 1)*eri%n_o + first_p,              &
                                   offset + (b-eri%n_o - 1)*eri%n_o + min(last_p, eri%n_o))
               enddo
            endif
         endif
!
         if (last_p .gt. eri%n_o) then !Read from vv
            offset = eri%n_o**2 + 2*eri%n_o*eri%n_v
            if(first_p .eq. eri%n_o+1 .and. last_p .eq. eri%n_mo) then
!
               call file_%read_(L_J_pq(:, :, max(first_q, eri%n_o+1):last_q), &
                                offset + (first_q - eri%n_o - 1)*eri%n_v + 1, &
                                offset + (last_q - eri%n_o)*eri%n_v)
            else
               do b = max(first_q, eri%n_o+1), last_q
!
                  call file_%read_(L_J_pq(:, max(first_p, eri%n_o+1):last_p, b),                 &
                                   offset + (b - eri%n_o-1)*eri%n_v + max(first_p - eri%n_o, 1), &
                                   offset + (b - eri%n_o-1)*eri%n_v + last_p - eri%n_o)
               enddo
            endif
         endif
      endif
!
   end subroutine read_cholesky_mo_eri_tool_c
!
!
   subroutine copy_in_cholesky_mo_eri_tool_c(eri, L_J_pq, &
                                             first_p, last_p, first_q, last_q, &
                                             L_J_oo, L_J_vo, L_J_ov, L_J_vv)
!!
!!    copy in cholesky
!!    Written by Rolf H. Myhre, May 2020
!!
!!    copy Cholesky vectors L_J_pq to arrays L_J_oo, L_J_vo, L_J_ov, and L_J_vv
!!
!!    first_p, last_p: first and last full index of p
!!    first_q, last_q: first and last full index of q
!!
!!    L_J_oo, L_J_vo, L_J_ov, L_J_vv: arrays to copy to
!!
      implicit none
!
      class(mo_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      complex(dp), dimension(eri%n_J, first_p:last_p, first_q:last_q), intent(in) :: L_J_pq
!
      complex(dp), dimension(eri%n_J, eri%n_o, eri%n_o), intent(out) :: L_J_oo
      complex(dp), dimension(eri%n_J, eri%n_v, eri%n_o), intent(out) :: L_J_vo
      complex(dp), dimension(eri%n_J, eri%n_o, eri%n_v), intent(out) :: L_J_ov
      complex(dp), dimension(eri%n_J, eri%n_v, eri%n_v), intent(out) :: L_J_vv
!
      integer :: i, j, a, b, K
!
!     oo block
!$omp parallel do private(j,i,K), collapse(2)
      do j = first_q, min(last_q, eri%n_o)
         do i = first_p, min(last_p, eri%n_o)
            do K = 1, eri%n_J
!
               L_J_oo(K, i, j) = L_J_pq(K, i, j)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     vo block
!$omp parallel do private(j,a,K), collapse(2)
      do j = first_q, min(last_q, eri%n_o)
         do a = max(first_p, eri%n_o+1), last_p
            do K = 1, eri%n_J
!
               L_J_vo(K, a - eri%n_o, j) = L_J_pq(K, a, j)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     ov block
!$omp parallel do private(b,i,K), collapse(2)
      do b = max(first_q, eri%n_o+1), last_q
         do i = first_p, min(last_p, eri%n_o)
            do K = 1, eri%n_J
!
               L_J_ov(K, i, b - eri%n_o) = L_J_pq(K, i, b)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     vv block
!$omp parallel do private(b,a,K), collapse(2)
      do b = max(first_q, eri%n_o+1), last_q
         do a = max(first_p, eri%n_o+1), last_p
            do K = 1, eri%n_J
!
               L_J_vv(K, a - eri%n_o, b - eri%n_o) = L_J_pq(K, a, b)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine copy_in_cholesky_mo_eri_tool_c
!
!
   subroutine copy_out_cholesky_mo_eri_tool_c(eri, L_J_pq, &
                                         first_p, last_p, first_q, last_q, &
                                         L_J_oo, L_J_vo, L_J_ov, L_J_vv)
!!
!!    copy out cholesky
!!    Written by Rolf H. Myhre, May 2020
!!
!!    copy Cholesky vectors from arrays L_J_oo, L_J_vo, L_J_ov, and L_J_vv to L_J_pq
!!
!!    first_p, last_p: first and last full index of p
!!    first_q, last_q: first and last full index of q
!!
!!    L_J_oo, L_J_vo, L_J_ov, L_J_vv: arrays to copy from
!!
      implicit none
!
      class(mo_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      complex(dp), dimension(eri%n_J, first_p:last_p, first_q:last_q), intent(out) :: L_J_pq
!
      complex(dp), dimension(eri%n_J, eri%n_o, eri%n_o), intent(in) :: L_J_oo
      complex(dp), dimension(eri%n_J, eri%n_v, eri%n_o), intent(in) :: L_J_vo
      complex(dp), dimension(eri%n_J, eri%n_o, eri%n_v), intent(in) :: L_J_ov
      complex(dp), dimension(eri%n_J, eri%n_v, eri%n_v), intent(in) :: L_J_vv
!
      integer :: i, j, a, b, K
!
!     oo block
!$omp parallel do private(j,i,K), collapse(2)
      do j = first_q, min(last_q, eri%n_o)
         do i = first_p, min(last_p, eri%n_o)
            do K = 1, eri%n_J
!
               L_J_pq(K, i, j) = L_J_oo(K, i, j)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     vo block
!$omp parallel do private(j,a,K), collapse(2)
      do j = first_q, min(last_q, eri%n_o)
         do a = max(first_p, eri%n_o+1), last_p
            do K = 1, eri%n_J
!
               L_J_pq(K, a, j) = L_J_vo(K, a - eri%n_o, j)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     ov block
!$omp parallel do private(b,i,K), collapse(2)
      do b = max(first_q, eri%n_o+1), last_q
         do i = first_p, min(last_p, eri%n_o)
            do K = 1, eri%n_J
!
               L_J_pq(K, i, b) = L_J_ov(K, i, b - eri%n_o)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     vv block
!$omp parallel do private(b,a,K), collapse(2)
      do b = max(first_q, eri%n_o+1), last_q
         do a = max(first_p, eri%n_o+1), last_p
            do K = 1, eri%n_J
!
               L_J_pq(K, a, b) = L_J_vv(K, a - eri%n_o, b - eri%n_o)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine copy_out_cholesky_mo_eri_tool_c
!
!
   subroutine set_cholesky_mo_mo_eri_tool_c(eri, L_J_pq, first_p, last_p, first_q, last_q)
!!
!!    Set Cholesky MO
!!    Written by Rolf H. Myhre, May 2020
!!
!!    based on routines by Sarai D. Folkestad and Eirik F. Kjønstad
!!
!!    Saves parts of the Cholesky MO vectors either to file or in memory.
!!
!!    L_J_pq:           Cholesky vector L_pq^J ordered as J x pq
!!
!!    first_p, last_p:  MO index range for p in L_J_pq
!!    first_q, last_q:  MO index range for q in L_J_pq
!!
      implicit none
!
      class(mo_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      complex(dp), dimension(eri%n_J, first_p:last_p, first_q:last_q), intent(in) :: L_J_pq
!
      if (last_p .lt. first_p) call output%error_msg('last_p .lt. first_p in set_cholesky_mo')
      if (last_q .lt. first_q) call output%error_msg('last_q .lt. first_q in set_cholesky_mo')
!
      if (eri%cholesky_mem) then
!
         call eri%copy_in_cholesky(L_J_pq, first_p, last_p, first_q, last_q, &
                                         eri%L_J_oo_mo, &
                                         eri%L_J_vo_mo, &
                                         eri%L_J_ov_mo, &
                                         eri%L_J_vv_mo)
      else
!
         call eri%cholesky_mo%open_('write')
!
         call eri%write_cholesky(L_J_pq, &
                                       first_p, last_p, &
                                       first_q, last_q, &
                                       eri%cholesky_mo)
!
         call eri%cholesky_mo%close_()
!
      endif
!
   end subroutine set_cholesky_mo_mo_eri_tool_c
!
!
   subroutine get_cholesky_mo_mo_eri_tool_c(eri, L_J_pq, first_p, last_p, first_q, last_q)
!!
!!    Get Cholesky MO
!!    Written by Rolf H. Myhre, May 2020
!!
!!    based on routines by Sarai D. Folkestad and Eirik F. Kjønstad
!!
!!    Returns parts of the Cholesky MO vectors either from file or memory.
!!
!!    L_J_pq:           Cholesky vector L_pq^J ordered as J x pq
!!
!!    first_p, last_p:  MO index range for p in L_J_pq
!!    first_q, last_q:  MO index range for q in L_J_pq
!!
      implicit none
!
      class(mo_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      complex(dp), dimension(eri%n_J, first_p:last_p, first_q:last_q), intent(out) :: L_J_pq
!
      if (last_p .lt. first_p) call output%error_msg('last_p .lt. first_p in get_cholesky_mo')
      if (last_q .lt. first_q) call output%error_msg('last_q .lt. first_q in get_cholesky_mo')
!
      if (eri%cholesky_mem) then
!
         call eri%copy_out_cholesky(L_J_pq, first_p, last_p, first_q, last_q, &
                                    eri%L_J_oo_mo, &
                                    eri%L_J_vo_mo, &
                                    eri%L_J_ov_mo, &
                                    eri%L_J_vv_mo)
      else
!
         call eri%cholesky_mo%open_('read')
!
         call eri%read_cholesky(L_J_pq, &
                                first_p, last_p, &
                                first_q, last_q, &
                                eri%cholesky_mo)
!
         call eri%cholesky_mo%close_()
!
      endif
!
   end subroutine get_cholesky_mo_mo_eri_tool_c
!
!
   subroutine update_cholesky_mo_mo_eri_tool_c(eri, T)
!!
!!    Update cholesky MO
!!    Written by Rolf H. Myhre, May 2020
!!
!!    based on routine written by Sarai D. Folkestad
!!
!!    Updates the Cholesky vectors on disk or memory
!!    using the transformation matrix T
!!
!!    L_J -> L'_J = T L_J T^T,
!!
!!    L''_J_rq = sum_s L_J_rs T_qs
!!    L'_J_pq = sum_r T_pr L''_J_rq
!!
!!    Temporary arrays L_J_1 and L_J_2 are used to store
!!    initial, intermediate and final vectors
!!
!!    Pointers L_J_1_p and L_J_2_p are used to keep track of changing dimensions.
!!
!!    In cases where we can't hold two full vectors in memory,
!!    we have to construct all the intermediate vectors in batches and write them to temp_file
!!    before reading them in and constructing the final vectors.
!!
      implicit none 
!
      class(mo_eri_tool_c), intent(inout) :: eri
!
      complex(dp), dimension(eri%n_mo, eri%n_mo), intent(in) :: T
!
      complex(dp), dimension(:), allocatable, target :: L_J_1,   L_J_2
      complex(dp), dimension(:,:,:), pointer         :: L_J_1_p, L_J_2_p
!
      type(batching_index) :: batcher
!
      logical :: all_in_mem
!
      integer :: batch, r
!
      type(direct_stream_file) :: temp_file
!
      batcher = batching_index(eri%n_mo)
!
      call mem%batch_setup(batcher, 0, 2*eri%n_J*eri%n_mo)
      all_in_mem = (batcher%num_batches .eq. 1)
!
      call mem%alloc(L_J_1, eri%n_J * batcher%max_length * eri%n_mo)
      call mem%alloc(L_J_2, eri%n_J * batcher%max_length * eri%n_mo)
!
      if (.not. all_in_mem) then
         temp_file = direct_stream_file('temp_file', eri%n_J, 2*dp, 'new')
         call temp_file%open_('readwrite')
      endif
!
!     Construct intermediate vectors L''_J_rq = sum_s L_J_rs T_qs
!
      do batch = 1, batcher%num_batches
!
         call batcher%determine_limits(batch)
!
         L_J_1_p(1:eri%n_J, 1:batcher%length, 1:eri%n_mo) => L_J_1
         L_J_2_p(1:eri%n_J, 1:batcher%length, 1:eri%n_mo) => L_J_2
!
         call eri%get_cholesky_mo(L_J_1, batcher%first, batcher%get_last(), 1, eri%n_mo)
!
         call zgemm('N', 'T',                                     &
                    eri%n_J*batcher%length, eri%n_mo, eri%n_mo,   &
                    one_complex,                                  &
                    L_J_1_p, eri%n_J*batcher%length,              &
                    T, eri%n_mo,                                  &
                    zero_complex,                                 &
                    L_J_2_p, eri%n_J*batcher%length)
!
         L_J_1_p(1:eri%n_J, 1:eri%n_mo, 1:batcher%length) => L_J_1
!
         call sort_123_to_132(L_J_2_p, L_J_1_p, eri%n_J, batcher%length, eri%n_mo)
!
         if (.not. all_in_mem) then
            call temp_file%write_(L_J_1, (batcher%first-1)*eri%n_mo+1, &
                                  batcher%get_last()*eri%n_mo)
         endif
!
      enddo
!
!     Construct final vectors L'_J_pq = sum_r T_pr L''_J_rq
!
      do batch = 1, batcher%num_batches
!
         call batcher%determine_limits(batch)
!
         L_J_1_p(1:eri%n_J, 1:batcher%length, 1:eri%n_mo) => L_J_1
         L_J_2_p(1:eri%n_J, 1:batcher%length, 1:eri%n_mo) => L_J_2
!
         if (.not. all_in_mem) then
            do r = 1, eri%n_mo
               call temp_file%read_(L_J_1_p(:, :, r), &
                                    (r-1)*eri%n_mo + batcher%first, &
                                    (r-1)*eri%n_mo + batcher%get_last())
            enddo
         endif
!
         call zgemm('N', 'T',                                     &
                    eri%n_J*batcher%length, eri%n_mo, eri%n_mo,   &
                    one_complex,                                  &
                    L_J_1_p, eri%n_J*batcher%length,              &
                    T, eri%n_mo,                                  &
                    zero_complex,                                 &
                    L_J_2_p, eri%n_J*batcher%length)

         L_J_1_p(1:eri%n_J, 1:eri%n_mo, 1:batcher%length) => L_J_1
!
         call sort_123_to_132(L_J_2_p, L_J_1_p, eri%n_J, batcher%length, eri%n_mo)
!
         call eri%set_cholesky_mo(L_J_1_p, 1, eri%n_mo, batcher%first, batcher%get_last())
!
      enddo

      if (.not. all_in_mem) then
         call temp_file%close_('delete')
      endif
!
      call mem%dealloc(L_J_1, eri%n_J * batcher%max_length * eri%n_mo)
      call mem%dealloc(L_J_2, eri%n_J * batcher%max_length * eri%n_mo)
!
   end subroutine update_cholesky_mo_mo_eri_tool_c
!
!
   pure subroutine get_eri_mo_mem_mo_eri_tool_c(eri, string, req_pq, req_rs, &
                                                dim_p, dim_q, dim_r, dim_s,  &
                                                qp, sr)
!!
!!    Get ERI MO mem
!!    Written by Rolf H. Myhre Jun 2020
!!
!!    Adds the memory usage for get ERI MO depending on pq and rs to req_pq and req_rs
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
      class(mo_eri_tool_c), intent(in) :: eri
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
!     No extra memory needed if integrals in mem
      if(.not. (eri%mo_eri_mem)) then
!
         full_p = ((string(1:1) .eq. 'v' .and. dim_p .eq. eri%n_v) .or. &
                   (string(1:1) .eq. 'o' .and. dim_p .eq. eri%n_o))
!
         full_r = ((string(3:3) .eq. 'v' .and. dim_r .eq. eri%n_v) .or. &
                   (string(3:3) .eq. 'o' .and. dim_r .eq. eri%n_o))
!
!        Can we use vectors in mem directly
         if(.not. (eri%cholesky_mem .and. full_p)) req_pq = req_pq + eri%n_J*dim_p*dim_q
         if(.not. (eri%cholesky_mem .and. full_r)) req_rs = req_rs + eri%n_J*dim_r*dim_s
!
!        Do we need temporary arrays for reordering
         if(present(qp)) then
            if(qp) req_pq = req_pq + eri%n_J*dim_p*dim_q
         endif
         if(present(sr)) then
            if(sr) req_rs = req_rs + eri%n_J*dim_r*dim_s
         endif
!
      endif
!
   end subroutine get_eri_mo_mem_mo_eri_tool_c
!
!
   subroutine get_eri_mo_mo_eri_tool_c(eri, string, g_pqrs, &
                                       first_p, last_p, first_q, last_q, &
                                       first_r, last_r, first_s, last_s, &
                                       alpha, beta, qp, sr, rspq)
!!
!!    Get ERI MO
!!    Written by Rolf H. Myhre Jun 2020
!!
!!    based on routines by Sarai D. Folkestad and Eirik F. Kjønstad
!!
!!    string: four character string indicating integral block. Example: "vovo"
!!            o: occupied, v: virtual, f: full
!!    q_pqrs: array to contain the integral
!!    first_p, last_p, etc.: first and last index of integrals in range determined by string
!!
!!    alpha: scales data added to g_pqrs (default = 1.0)
!!    beta : scales data already in g_pqrs (default = 0.0)
!!
!!    Optional reordering logicals (default = .false.):
!!    qp:   switches order of p and q, i.e., g_qprs
!!    sr:   switches order of r and s, i.e., g_pqsr
!!    rspq: switches pq and rs, i.e., g_rspq
!!    All these can be combined freely for eight different combinations
!!
      implicit none
!
      class(mo_eri_tool_c), intent(inout) :: eri
!
      character(len=4), intent(in) :: string
!
      integer, optional, intent(in) :: first_p, last_p
      integer, optional, intent(in) :: first_q, last_q
      integer, optional, intent(in) :: first_r, last_r
      integer, optional, intent(in) :: first_s, last_s
!
      complex(dp), intent(inout), dimension(1) :: g_pqrs
!
      complex(dp), optional, intent(in) :: alpha, beta
!
      logical, optional, intent(in) :: qp, sr, rspq
!
      integer :: full_first_p, full_last_p, full_first_q, full_last_q
      integer :: full_first_r, full_last_r, full_first_s, full_last_s
!
      complex(dp) :: alpha_, beta_
!
      alpha_ = one_complex
      if(present(alpha)) alpha_ = alpha
!
      beta_ = zero_complex
      if(present(beta)) beta_ = beta
!
      call eri%index_setup(string(1:1), full_first_p, full_last_p, first_p, last_p)
      call eri%index_setup(string(2:2), full_first_q, full_last_q, first_q, last_q)
      call eri%index_setup(string(3:3), full_first_r, full_last_r, first_r, last_r)
      call eri%index_setup(string(4:4), full_first_s, full_last_s, first_s, last_s)
!
      call eri%get_g_pqrs_mo(g_pqrs, full_first_p, full_last_p, &
                                     full_first_q, full_last_q, &
                                     full_first_r, full_last_r, &
                                     full_first_s, full_last_s, &
                                     alpha_, beta_, qp, sr, rspq)
!
   end subroutine get_eri_mo_mo_eri_tool_c
!
!
   subroutine get_g_pqrs_mo_mo_eri_tool_c(eri, g_pqrs, first_p, last_p, first_q, last_q, &
                                                       first_r, last_r, first_s, last_s, &
                                                       alpha, beta, qp, sr, rspq)
!!
!!    Get g_pqrs MO
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    based on routines by Sarai D. Folkestad and Eirik F. Kjønstad
!!
!!    q_pqrs: array to contain the integral
!!    first_p, last_p, etc.: absolute first and last index of integrals
!!
!!    Optional reordering logicals (default = .false.):
!!    qp:   switches order of p and q, i.e., g_qprs
!!    sr:   switches order of r and s, i.e., g_pqsr
!!    rspq: switches pq and rs, i.e., g_rspq
!!    All these can be combined freely for eight different combinations
!!
      implicit none
!
      class(mo_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
      integer, intent(in) :: first_r, last_r
      integer, intent(in) :: first_s, last_s
!
      complex(dp), intent(inout), &
                dimension(first_p:last_p,first_q:last_q,first_r:last_r,first_s:last_s) :: g_pqrs
!
      complex(dp), intent(in) :: alpha, beta
!
      logical, optional, intent(in) :: qp, sr, rspq
!
      if (eri%mo_eri_mem) then
!
         call eri%copy_g_pqrs(g_pqrs, eri%g_pqrs_mo, &
                              first_p, last_p, first_q, last_q, &
                              first_r, last_r, first_s, last_s, &
                              alpha, beta, qp, sr, rspq)
!
      else
!
         call eri%construct_g_pqrs_mo(g_pqrs, first_p, last_p, first_q, last_q, &
                                              first_r, last_r, first_s, last_s, &
                                              alpha, beta, qp, sr, rspq)
!
      endif
!
   end subroutine get_g_pqrs_mo_mo_eri_tool_c
!
!
   subroutine construct_g_pqrs_mo_mo_eri_tool_c(eri, g_pqrs, &
                                           first_p, last_p, first_q, last_q, &
                                           first_r, last_r, first_s, last_s, &
                                           alpha, beta, qp, sr, rspq)
!!
!!    Construct g_pqrs MO
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    q_pqrs: array to contain the integral
!!    first_p, last_p, etc.: absolute first and last index of integrals
!!
!!    Optional reordering logicals (default = .false.):
!!    qp:   switches order of p and q, i.e., g_qprs
!!    sr:   switches order of r and s, i.e., g_pqsr
!!    rspq: switches pq and rs, i.e., g_rspq
!!    All these can be combined freely for eight different combinations
!!
      implicit none
!
      class(mo_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
      integer, intent(in) :: first_r, last_r
      integer, intent(in) :: first_s, last_s
!
      complex(dp), intent(inout), &
                dimension(first_p:last_p,first_q:last_q,first_r:last_r,first_s:last_s) :: g_pqrs
!
      complex(dp), intent(in) :: alpha, beta
!
      logical, optional, intent(in) :: qp, sr, rspq
!
      complex(dp), dimension(:,:,:), pointer, contiguous :: L_J_pq_p, L_J_rs_p
!
      integer :: dim_p, dim_q, dim_r, dim_s
!
      logical :: switch_pq, switch_rs
!
      logical :: pq_alloced, rs_alloced
!
      L_J_pq_p => null()
      L_J_rs_p => null()
!
      switch_pq = .false.
      switch_rs = .false.
      if(present(qp)) switch_pq = qp
      if(present(sr)) switch_rs = sr
!
      dim_p = last_p - first_p + 1
      dim_q = last_q - first_q + 1
      dim_r = last_r - first_r + 1
      dim_s = last_s - first_s + 1
!
      call eri%L_pointer_setup_mo(L_J_pq_p, switch_pq, &
                                  first_p, last_p, first_q, last_q, pq_alloced)
!
      if((first_p .eq. first_r) .and. (last_p .eq. last_r) .and. &
         (first_q .eq. first_s) .and. (last_q .eq. last_s) .and. &
         (switch_pq .eqv. switch_rs)) then
!
         call eri%construct_g_symm_from_L(L_J_pq_p, g_pqrs, alpha, beta, dim_p, dim_q)
!
      else
         call eri%L_pointer_setup_mo(L_J_rs_p, switch_rs, &
                                     first_r, last_r, first_s, last_s, rs_alloced)
!
         call eri%construct_g_from_L(L_J_pq_p, L_J_rs_p, g_pqrs, alpha, beta, &
                                     dim_p, dim_q, dim_r, dim_s, rspq)
!
         if (rs_alloced) then
            call mem%dealloc(L_J_rs_p, eri%n_J, dim_r, dim_s)
         endif
      endif
!
      if (pq_alloced) then
         call mem%dealloc(L_J_pq_p, eri%n_J, dim_p, dim_q)
      endif
!
   end subroutine construct_g_pqrs_mo_mo_eri_tool_c
!
!
   subroutine copy_g_pqrs_mo_eri_tool_c(eri, g_to, g_from, &
                                        first_p, last_p, first_q, last_q, &
                                        first_r, last_r, first_s, last_s, &
                                        alpha, beta, qp, sr, rspq)
!!
!!    Copy g_pqrs MO
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    g_to: array to copy integrals to
!!    g_from: array to copy integrals from
!!    first_p, last_p, etc.: absolute first and last index of integrals
!!
!!    Optional reordering logicals (default = .false.):
!!    qp:   switches order of p and q, i.e., g_qprs
!!    sr:   switches order of r and s, i.e., g_pqsr
!!    rspq: switches pq and rs, i.e., g_rspq
!!    All these can be combined freely for eight different combinations
!!
      implicit none
!
      class(mo_eri_tool_c), intent(in) :: eri
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
      integer, intent(in) :: first_r, last_r
      integer, intent(in) :: first_s, last_s
!
      complex(dp), target, intent(inout), &
             dimension(first_p:last_p, first_q:last_q, first_r:last_r, first_s:last_s) :: g_to
!
      complex(dp), intent(in), dimension(eri%n_mo, eri%n_mo, eri%n_mo, eri%n_mo) :: g_from
!
      complex(dp), intent(in) :: alpha, beta
!
      logical, optional, intent(in) :: qp, sr, rspq
!
      integer :: p, q, r, s
      logical :: qp_, sr_, rspq_
!
!     We need a pointer to get the right dimensions in case indices are switched
!
      complex(dp), dimension(:,:,:,:), pointer :: g_to_p => null()
!
      qp_ = .false.
      sr_ = .false.
      rspq_ = .false.
      if(present(qp)) qp_ = qp
      if(present(sr)) sr_ = sr
      if(present(rspq)) rspq_ = rspq
!
      if (beta .eq. zero_complex) then
!
         if(.not. qp_ .and. .not. sr_ .and. .not. rspq_) then !pqrs
!
!$omp parallel do private(s, r, q, p)
            do s = first_s, last_s
               do r = first_r, last_r
                  do q = first_q, last_q
                     do p = first_p, last_p
                        g_to(p,q,r,s) = alpha*g_from(p,q,r,s)
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         elseif(qp_ .and. .not. sr_ .and. .not. rspq_) then !qprs
!
            g_to_p(first_q:last_q, first_p:last_p, first_r:last_r, first_s:last_s) => g_to
!
!$omp parallel do private(s, r, q, p)
            do s = first_s, last_s
               do r = first_r, last_r
                  do p = first_p, last_p
                     do q = first_q, last_q
                        g_to_p(q,p,r,s) = alpha*g_from(p,q,r,s)
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         elseif(.not. qp_ .and. sr_ .and. .not. rspq_) then !pqsr
!
            g_to_p(first_p:last_p, first_q:last_q, first_s:last_s, first_r:last_r) => g_to
!
!$omp parallel do private(s, r, q, p)
            do r = first_r, last_r
               do s = first_s, last_s
                  do q = first_q, last_q
                     do p = first_p, last_p
                        g_to_p(p,q,s,r) = alpha*g_from(p,q,r,s)
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         elseif(qp_ .and. sr_ .and. .not. rspq_) then !qpsr
!
            g_to_p(first_q:last_q, first_p:last_p, first_s:last_s, first_r:last_r) => g_to
!
!$omp parallel do private(s, r, q, p)
            do r = first_r, last_r
               do s = first_s, last_s
                  do p = first_p, last_p
                     do q = first_q, last_q
                        g_to_p(q,p,s,r) = alpha*g_from(p,q,r,s)
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         elseif(.not. qp_ .and. .not. sr_ .and. rspq_) then !rspq
!
            g_to_p(first_r:last_r, first_s:last_s, first_p:last_p, first_q:last_q) => g_to
!
!$omp parallel do private(s, r, q, p)
            do q = first_q, last_q
               do p = first_p, last_p
                  do s = first_s, last_s
                     do r = first_r, last_r
                        g_to_p(r,s,p,q) = alpha*g_from(p,q,r,s)
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         elseif(qp_ .and. .not. sr_ .and. rspq_) then !rsqp
!
            g_to_p(first_r:last_r, first_s:last_s, first_q:last_q, first_p:last_p) => g_to
!
!$omp parallel do private(s, r, q, p)
            do p = first_p, last_p
               do q = first_q, last_q
                  do s = first_s, last_s
                     do r = first_r, last_r
                        g_to_p(r,s,q,p) = alpha*g_from(p,q,r,s) !srpq
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         elseif(.not. qp_ .and. sr_ .and. rspq_) then !srpq
!
            g_to_p(first_s:last_s, first_r:last_r, first_p:last_p, first_q:last_q) => g_to
!
!$omp parallel do private(s, r, q, p)
            do q = first_q, last_q
               do p = first_p, last_p
                  do r = first_r, last_r
                     do s = first_s, last_s
                        g_to_p(s,r,p,q) = alpha*g_from(p,q,r,s)
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         elseif(qp_ .and. sr_ .and. rspq_) then !srqp
!
            g_to_p(first_s:last_s, first_r:last_r, first_q:last_q, first_p:last_p) => g_to
!
!$omp parallel do private(s, r, q, p)
            do p = first_p, last_p
               do q = first_q, last_q
                  do r = first_r, last_r
                     do s = first_s, last_s
                        g_to_p(s,r,q,p) = alpha*g_from(p,q,r,s)
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         endif
!
      else !beta .ne. zero
!
         if(.not. qp_ .and. .not. sr_ .and. .not. rspq_) then !pqrs
!
!$omp parallel do private(s, r, q, p)
            do s = first_s, last_s
               do r = first_r, last_r
                  do q = first_q, last_q
                     do p = first_p, last_p
                        g_to(p,q,r,s) = beta*g_to(p,q,r,s) + alpha*g_from(p,q,r,s)
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         elseif(qp_ .and. .not. sr_ .and. .not. rspq_) then !qprs
!
            g_to_p(first_q:last_q, first_p:last_p, first_r:last_r, first_s:last_s) => g_to
!
!$omp parallel do private(s, r, q, p)
            do s = first_s, last_s
               do r = first_r, last_r
                  do p = first_p, last_p
                     do q = first_q, last_q
                        g_to_p(q,p,r,s) = beta*g_to_p(q,p,r,s) + alpha*g_from(p,q,r,s)
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         elseif(.not. qp_ .and. sr_ .and. .not. rspq_) then !pqsr
!
            g_to_p(first_p:last_p, first_q:last_q, first_s:last_s, first_r:last_r) => g_to
!
!$omp parallel do private(s, r, q, p)
            do r = first_r, last_r
               do s = first_s, last_s
                  do q = first_q, last_q
                     do p = first_p, last_p
                        g_to_p(p,q,s,r) = beta*g_to_p(p,q,s,r) + alpha*g_from(p,q,r,s)
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         elseif(qp_ .and. sr_ .and. .not. rspq_) then !qpsr
!
            g_to_p(first_q:last_q, first_p:last_p, first_s:last_s, first_r:last_r) => g_to
!
!$omp parallel do private(s, r, q, p)
            do r = first_r, last_r
               do s = first_s, last_s
                  do p = first_p, last_p
                     do q = first_q, last_q
                        g_to_p(q,p,s,r) = beta*g_to_p(q,p,s,r) + alpha*g_from(p,q,r,s)
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         elseif(.not. qp_ .and. .not. sr_ .and. rspq_) then !rspq
!
            g_to_p(first_r:last_r, first_s:last_s, first_p:last_p, first_q:last_q) => g_to
!
!$omp parallel do private(s, r, q, p)
            do q = first_q, last_q
               do p = first_p, last_p
                  do s = first_s, last_s
                     do r = first_r, last_r
                        g_to_p(r,s,p,q) = beta*g_to_p(r,s,p,q) + alpha*g_from(p,q,r,s)
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         elseif(qp_ .and. .not. sr_ .and. rspq_) then !rsqp
!
            g_to_p(first_r:last_r, first_s:last_s, first_q:last_q, first_p:last_p) => g_to
!
!$omp parallel do private(s, r, q, p)
            do p = first_p, last_p
               do q = first_q, last_q
                  do s = first_s, last_s
                     do r = first_r, last_r
                        g_to_p(r,s,q,p) = beta*g_to_p(r,s,q,p) + alpha*g_from(p,q,r,s) !srpq
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         elseif(.not. qp_ .and. sr_ .and. rspq_) then !srpq
!
            g_to_p(first_s:last_s, first_r:last_r, first_p:last_p, first_q:last_q) => g_to
!
!$omp parallel do private(s, r, q, p)
            do q = first_q, last_q
               do p = first_p, last_p
                  do r = first_r, last_r
                     do s = first_s, last_s
                        g_to_p(s,r,p,q) = beta*g_to_p(s,r,p,q) + alpha*g_from(p,q,r,s)
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         elseif(qp_ .and. sr_ .and. rspq_) then !srqp
!
            g_to_p(first_s:last_s, first_r:last_r, first_q:last_q, first_p:last_p) => g_to
!
!$omp parallel do private(s, r, q, p)
            do p = first_p, last_p
               do q = first_q, last_q
                  do r = first_r, last_r
                     do s = first_s, last_s
                        g_to_p(s,r,q,p) = beta*g_to_p(s,r,q,p) + alpha*g_from(p,q,r,s)
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         endif
!
      endif !beta .eq. zero
!
   end subroutine copy_g_pqrs_mo_eri_tool_c
!
!
   subroutine copy_g_pqrs_to_packed_mo_eri_tool_c(eri, g_to, g_from, &
                                               first_p, first_q, &
                                               dim_p, dim_q, alpha, beta, qp)
!!
!!    Copy packed g_pqrs
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    g_to: 2d array in RFP format to copy integrals to
!!    g_from: 4d array to copy integrals from
!!    first_p, first_q.: absolute first indices of integrals
!!    dim_p, dim_q: length of p and q range
!!
!!    Optional reordering logicals (default = .false.):
!!    qp:   switches order of p and q, i.e., reorders to g_qpqp
!!
      implicit none
!
      class(mo_eri_tool_c), intent(in) :: eri
!
      integer, intent(in) :: first_p, first_q
      integer, intent(in) :: dim_p, dim_q
!
      complex(dp), intent(inout), &
                   dimension(dim_p*dim_q + mod(dim_p*dim_q+1,2), (dim_p*dim_q+1)/2) :: g_to
!
      complex(dp), intent(in), dimension(eri%n_mo, eri%n_mo, eri%n_mo, eri%n_mo) :: g_from
!
      complex(dp), intent(in) :: alpha, beta
!
      logical, optional, intent(in) :: qp
!
      integer :: p, q, r, s, pq, rs, x, y
      integer :: tridim
      logical :: qp_
!
      tridim = dim_p*dim_q/2
!
      qp_ = .false.
      if(present(qp)) qp_ = qp
!
      if(beta .eq. zero_complex) then
!
         if(.not. qp_) then !pqpq
!
!$omp parallel do private(x, y, pq, rs, p, q, r, s)
            do y = 1, (dim_p*dim_q+1)/2
               do x = 1, dim_p*dim_q + mod(dim_p*dim_q+1,2)
!         
                  if (x .le. y + tridim) then
                     pq = x
                     rs = y + tridim
                  else
                     pq = y
                     rs = x - tridim - 1
                  endif
!         
                  p = mod(pq-1,dim_p)
                  r = mod(rs-1,dim_p)
                  q = (pq-1)/dim_p
                  s = (rs-1)/dim_p
!         
                  g_to(x,y) = alpha*g_from(first_p+p, first_q+q, first_p+r, first_q+s)
!         
               enddo
            enddo
!$omp end parallel do
!
         else !qpqp
!
!$omp parallel do private(x, y, pq, rs, p, q, r, s)
            do y = 1, (dim_p*dim_q+1)/2
               do x = 1, dim_p*dim_q + mod(dim_p*dim_q+1,2)
!        
                  if (x .le. y + tridim) then
                     pq = x
                     rs = y + tridim
                  else
                     pq = y
                     rs = x - tridim - 1
                  endif
!        
                  q = mod(pq-1,dim_q)
                  s = mod(rs-1,dim_q)
                  p = (pq-1)/dim_q
                  r = (rs-1)/dim_q
!        
                  g_to(x,y) = alpha*g_from(first_p+p, first_q+q, first_p+r, first_q+s)
!        
               enddo
            enddo
!$omp end parallel do
!
         endif
!
      else !beta .ne. zero
!
         if(.not. qp_) then !pqpq
!
!$omp parallel do private(x, y, pq, rs, p, q, r, s)
            do y = 1, (dim_p*dim_q+1)/2
               do x = 1, dim_p*dim_q + mod(dim_p*dim_q+1,2)
!         
                  if (x .le. y + tridim) then
                     pq = x
                     rs = y + tridim
                  else
                     pq = y
                     rs = x - tridim - 1
                  endif
!         
                  p = mod(pq-1,dim_p)
                  r = mod(rs-1,dim_p)
                  q = (pq-1)/dim_p
                  s = (rs-1)/dim_p
!         
                  g_to(x,y) = beta*g_to(x,y) &
                            + alpha*g_from(first_p+p, first_q+q, first_p+r, first_q+s)
!         
               enddo
            enddo
!$omp end parallel do
!
         else !qpqp
!
!$omp parallel do private(x, y, pq, rs, p, q, r, s)
            do y = 1, (dim_p*dim_q+1)/2
               do x = 1, dim_p*dim_q + mod(dim_p*dim_q+1,2)
!        
                  if (x .le. y + tridim) then
                     pq = x
                     rs = y + tridim
                  else
                     pq = y
                     rs = x - tridim - 1
                  endif
!        
                  q = mod(pq-1,dim_q)
                  s = mod(rs-1,dim_q)
                  p = (pq-1)/dim_q
                  r = (rs-1)/dim_q
!        
                  g_to(x,y) = beta*g_to(x,y) &
                            + alpha*g_from(first_p+p, first_q+q, first_p+r, first_q+s)
!        
               enddo
            enddo
!$omp end parallel do
!
         endif
      endif !beta .eq. zero
!
   end subroutine copy_g_pqrs_to_packed_mo_eri_tool_c
!
!
   subroutine place_g_mo_in_memory_mo_eri_tool_c(eri)
!!
!!    Place g MO in memory
!!    Written by Rolf H. Myhre, Aug 2020
!!
      implicit none
!
      class(mo_eri_tool_c), intent(inout) :: eri
!
      if (eri%mo_eri_mem) then
!
         call eri%construct_g_pqrs_mo(eri%g_pqrs_mo, 1, eri%n_mo, 1, eri%n_mo, &
                                                     1, eri%n_mo, 1, eri%n_mo, &
                                                     one_complex, zero_complex)
!
      endif
!
   end subroutine place_g_mo_in_memory_mo_eri_tool_c
!
!
   subroutine construct_g_from_L_mo_eri_tool_c(eri, L_J_pq, L_J_rs, g_pqrs, alpha, beta, &
                                               dim_p, dim_q, dim_r, dim_s, rspq)
!!
!!    construct g from L
!!    written by Rolf H. Myhre, Jun 2020
!!
!!    Contract the Cholesky vectors L_J_pq and L_J_rs to construct g_pqrs or g_rspq
!!
      implicit none
!
      class(mo_eri_tool_c), intent(in) :: eri
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(eri%n_J*dim_p*dim_q), intent(in) :: L_J_pq
      complex(dp), dimension(eri%n_J*dim_r*dim_s), intent(in) :: L_J_rs
!
      complex(dp), dimension(dim_p*dim_q*dim_r*dim_s), intent(inout) :: g_pqrs
!
      complex(dp), intent(in) :: alpha, beta
!
      logical, optional, intent(in) :: rspq
      logical :: switch_pq_rs
!
      switch_pq_rs = .false.
      if (present(rspq)) switch_pq_rs = rspq
!
      if (.not. switch_pq_rs) then
         call zgemm('T', 'N',     &
                     dim_p*dim_q, &
                     dim_r*dim_s, &
                     eri%n_J,     &
                     alpha,       &
                     L_J_pq,      &
                     eri%n_J,     &
                     L_J_rs,      &
                     eri%n_J,     &
                     beta,        &
                     g_pqrs,      &
                     dim_p*dim_q)
      else
         call zgemm('T', 'N',     &
                     dim_r*dim_s, &
                     dim_p*dim_q, &
                     eri%n_J,     &
                     alpha,       &
                     L_J_rs,      &
                     eri%n_J,     &
                     L_J_pq,      &
                     eri%n_J,     &
                     beta,        &
                     g_pqrs,      &
                     dim_r*dim_s)
      endif
!
   end subroutine construct_g_from_L_mo_eri_tool_c
!
!
   subroutine construct_g_symm_from_L_mo_eri_tool_c(eri, L_J_pq, g_pqpq, alpha, beta, &
                                                    dim_p, dim_q)
!!
!!    construct packed g from L
!!    written by Rolf H. Myhre, Jan 2021
!!
!!    Contract the Cholesky vector L_J_pq to construct g_pqpq
!!
      implicit none
!
      class(mo_eri_tool_c), intent(in) :: eri
!
      integer, intent(in) :: dim_p, dim_q
!
      complex(dp), dimension(eri%n_J, dim_p*dim_q), intent(in) :: L_J_pq
!
      complex(dp), intent(inout), dimension(dim_p*dim_q, dim_p*dim_q) :: g_pqpq
!
      complex(dp), intent(in) :: alpha, beta
!
      integer :: p, q
!
      call zsyrk('U', 'T', dim_p*dim_q, eri%n_J, alpha, L_J_pq, eri%n_J, beta, g_pqpq, dim_p*dim_q)
!
      if(beta .eq. zero) then
!$omp parallel do private(p,q)
         do q = 1, dim_p*dim_q
            do p = q+1, dim_p*dim_q
               g_pqpq(p,q) = g_pqpq(q,p)
            enddo
         enddo
!$omp end parallel do
      else
!$omp parallel do private(p,q)
         do q = 1, dim_p*dim_q
            do p = q+1, dim_p*dim_q
               g_pqpq(p,q) = beta*g_pqpq(p,q) + g_pqpq(q,p)
            enddo
         enddo
!$omp end parallel do
      endif
!
   end subroutine construct_g_symm_from_L_mo_eri_tool_c
!
!
   subroutine construct_g_packed_from_L_mo_eri_tool_c(eri, L_J_pq, g_pqpq, alpha, beta, &
                                                      dim_p, dim_q)
!!
!!    construct packed g from L
!!    written by Rolf H. Myhre, Jan 2021
!!
!!    Contract the Cholesky vectors L_J_pq to construct g_pqpq
!!
!!    Note, no zsrfk routine, so have to call zsyrk twice and dgemm once
!!
      implicit none
!
      class(mo_eri_tool_c), intent(in) :: eri
!
      integer, intent(in) :: dim_p, dim_q
!
      complex(dp), dimension(eri%n_J, dim_p*dim_q), intent(in) :: L_J_pq
!
      complex(dp), intent(inout), dimension(dim_p*dim_q*(dim_p*dim_q+1)/2) :: g_pqpq
!
      complex(dp), intent(in) :: alpha, beta
!
      integer :: n, n1, n2, off
!
      n = dim_p*dim_q
      n1 = n/2
      n2 = n - n1
      off = mod(n+1,2)
!
      call zsyrk('L', 'T', n1, eri%n_J, alpha, L_J_pq(1, 1), eri%n_J, &
                 beta, g_pqpq(n2+1+off), n+off)
!
      call zsyrk('U', 'T', n2, eri%n_J, alpha, L_J_pq(1, n2+off), eri%n_J, &
                 beta, g_pqpq(n1+1), n+off)
!
      call zgemm('T', 'N', n1, n2, eri%n_J, alpha, L_J_pq(1, 1), &
                 eri%n_J, L_J_pq(1, n2+off), eri%n_J, beta, g_pqpq(1), n+off)
!
   end subroutine construct_g_packed_from_L_mo_eri_tool_c
!
!
   subroutine set_pointer_mo_mo_eri_tool_c(eri, point, last_p, last_q, first_q)
!!
!!    set pointer MO
!!    written by Rolf H. Myhre, Jun 2020
!!
!!    Set pointer to Cholesky MO vector block
!!    Assumes full index p
!!
      implicit none
!
      class(mo_eri_tool_c), intent(in), target :: eri
      complex(dp), dimension(:,:,:), pointer :: point
      integer, intent(in) :: last_p, last_q, first_q
!
      if(associated(point)) call output%error_msg('Pointer associated in set pointer mo')
      if(.not. eri%cholesky_mem) call output%error_msg('Vectors not in mem for set poiner mo')
!
      if (last_q .le. eri%n_o) then
         if (last_p .le. eri%n_o) then
            point => eri%L_J_oo_mo(:,:,first_q:last_q)
         else
            point => eri%L_J_vo_mo(:,:,first_q:last_q)
         endif
      else
         if (last_p .le. eri%n_o) then
            point => eri%L_J_ov_mo(:,:,first_q-eri%n_o:last_q-eri%n_o)
         else
            point => eri%L_J_vv_mo(:,:,first_q-eri%n_o:last_q-eri%n_o)
         endif
      endif
!
   end subroutine set_pointer_mo_mo_eri_tool_c
!
!
   subroutine L_pointer_setup_mo_mo_eri_tool_c(eri, point, switch, &
                                               first_p, last_p, first_q, last_q, alloced)
!!
!!    pointer setup MO
!!    written by Rolf H. Myhre, Jun 2020
!!
!!    Do not use this routine!
!!
!!    Sets the pointer 'point' to either point to a stored Cholesky vector
!!    or allocates the pointer and fills it with a vector in appropriate order
!!
!!    switch: logical, 'point' will point to reordered vectors if .true.
!!            that is, L is J,q,p ordered instead of J,p,q ordered
!!
!!    Deallocing 'point' if alloced = .false. will take you on exciting adventures
!!
      implicit none
!
      class(mo_eri_tool_c), intent(inout) :: eri
      complex(dp), dimension(:,:,:), pointer, contiguous, intent(out) :: point
      logical, intent(in) :: switch
      integer, intent(in) :: first_p, last_p, first_q, last_q
      logical, intent(out) :: alloced
!
      complex(dp), dimension(:,:,:), pointer, contiguous :: temp
      logical :: vec_in_mem
!
      if(associated(point)) call output%error_msg('Pointer associated in pointer setup')

      temp => null()
!
!     Check if vector is in memory and contiguous in a single block
!
      vec_in_mem = (eri%cholesky_mem .and. &
                    eri%is_pq_contiguous(first_p, last_p, first_q, last_q))
!
!     Allocate point if vector is not in memory and/or we are reordering
!
      alloced = .false.
      if(.not. vec_in_mem) then
         call mem%alloc(point, eri%n_J, last_p-first_p+1, last_q-first_q+1)
         alloced = .true.
      elseif(switch) then
         call mem%alloc(point, eri%n_J, last_q-first_q+1, last_p-first_p+1)
         alloced = .true.
      endif
!
!     Allocate temporary array if vector not in memory and reordering
!
      if(.not. vec_in_mem .and. switch) then
         call mem%alloc(temp, eri%n_J, last_p-first_p+1, last_q-first_q+1)
      endif
!
      if(.not. switch) then !No Reordering
!
!        Set point directly if possible, else read from file
!
         if(vec_in_mem) then
            call eri%set_pointer_mo(point, last_p, last_q, first_q)
         else
            call eri%get_cholesky_mo(point, first_p, last_p, first_q, last_q)
         endif
!
      else !Reordering
!
!        Set temp to vectors if possible, else read from file
!
         if(vec_in_mem) then
            call eri%set_pointer_mo(temp, last_p, last_q, first_q)
         else
            call eri%get_cholesky_mo(temp, first_p, last_p, first_q, last_q)
         endif
!
!        Sort into point
!
         call sort_123_to_132(temp, point, eri%n_J, last_p-first_p+1, last_q-first_q+1)
!
!        Unset temp if in memory, else dealloc
!
         if (vec_in_mem) then
            temp => null()
         else
            call mem%dealloc(temp, eri%n_J, last_p-first_p+1, last_q-first_q+1)
         endif
!
      endif
!
   end subroutine L_pointer_setup_mo_mo_eri_tool_c
!
!
   subroutine cleanup_mo_eri_tool_c(eri)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2019
!!
      implicit none
!
      class(mo_eri_tool_c) :: eri
!
      if (eri%cholesky_mem) then
!
!        Deallocate Cholesky arrays
!
         call mem%dealloc(eri%L_J_oo_mo, eri%n_J, eri%n_o, eri%n_o)
         call mem%dealloc(eri%L_J_ov_mo, eri%n_J, eri%n_o, eri%n_v)
         call mem%dealloc(eri%L_J_vo_mo, eri%n_J, eri%n_v, eri%n_o)
         call mem%dealloc(eri%L_J_vv_mo, eri%n_J, eri%n_v, eri%n_v)
!
      endif
!
      if (eri%mo_eri_mem) then
!
!        Deallocate electron repulsion integrals
!
         call mem%dealloc(eri%g_pqrs_mo, eri%n_mo, eri%n_mo, &
                                         eri%n_mo, eri%n_mo)
!
      endif
!
      eri%cholesky_mem = .false.
      eri%mo_eri_mem = .false.
!
   end subroutine cleanup_mo_eri_tool_c
!
!
end module mo_eri_tool_c_class
