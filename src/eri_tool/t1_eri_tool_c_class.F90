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
module t1_eri_tool_c_class
!!
!!    T1 ERI tool complex class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Modified by Rolf H. Myhre, June 2020
!!
!!    Tool that handles the two-electron repulsion integrals in the MO/T1/"C1"
!!    basis. The tool uses Cholesky vectors L_pq^J to construct
!!    the electron repulsion integrals:
!!
!!       g_pqrs = sum_J L_pq^J L_rs^J
!!
!!    set_t1_to_mo copies the MO Cholesky vectors to the T1 vectors
!!    update_t1_integrals constructs the T1 Cholesky vectors using the provided T1 amplitudes
!!
!!    Calling get_eri_t1 returns the integrals
!!    Calling get_eri_t1_mem returns an estimate for the amount of memory needed
!!
!
   use global_out,   only : output
   use global_in,    only : input
!
   use mo_eri_tool_class,        only : mo_eri_tool
   use t1_eri_tool_class,        only : t1_eri_tool
   use mo_eri_tool_c_class,      only : mo_eri_tool_c
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
   type, extends(mo_eri_tool_c) :: t1_eri_tool_c
!
!     Files to keep L_pq^J in MO and T1 bases
!
      type(direct_stream_file) :: cholesky_t1
!
!     Arrays for T1 Cholesky vectors
!
      complex(dp), dimension(:,:,:), allocatable :: L_J_oo_t1
      complex(dp), dimension(:,:,:), allocatable :: L_J_vo_t1
      complex(dp), dimension(:,:,:), allocatable :: L_J_ov_t1
      complex(dp), dimension(:,:,:), allocatable :: L_J_vv_t1
!
      logical :: t1_eri_mem
      complex(dp), dimension(:,:,:,:), allocatable :: g_pqrs_t1
!
   contains
!
!     Read, print settings
!
      procedure :: print_settings            => print_settings_t1_eri_tool_c
      procedure :: read_settings             => read_settings_t1_eri_tool_c
!
!     Cleanup
!
      procedure :: cleanup                   => cleanup_t1_eri_tool_c
!
!     Initialization
!
      procedure :: initialize                => initialize_t1_eri_tool_c
      procedure :: copy_from_mo              => copy_from_mo_t1_eri_tool_c
      procedure :: copy_from_real_mo         => copy_from_real_mo_t1_eri_tool_c
      procedure :: copy_from_t1              => copy_from_t1_t1_eri_tool_c
      procedure :: copy_from_real_t1         => copy_from_real_t1_t1_eri_tool_c
!
!     T1 Cholesky vectors
!
      procedure :: get_cholesky_t1           => get_cholesky_t1_t1_eri_tool_c
      procedure :: set_cholesky_t1           => set_cholesky_t1_t1_eri_tool_c
!
      procedure :: set_t1_to_mo              => set_t1_to_mo_t1_eri_tool_c
      procedure :: update_t1_integrals       => update_t1_integrals_t1_eri_tool_c
!
      procedure :: construct_t1_cholesky     => construct_t1_cholesky_t1_eri_tool_c
      procedure :: construct_cholesky_t1_oo  => construct_cholesky_t1_oo_t1_eri_tool_c
      procedure :: construct_cholesky_t1_vo  => construct_cholesky_t1_vo_t1_eri_tool_c
      procedure :: construct_cholesky_t1_vv  => construct_cholesky_t1_vv_t1_eri_tool_c
!
!     Get and construct routines for g_pqrs in T1 basis
!
      procedure :: get_eri_t1_mem            => get_eri_t1_mem_t1_eri_tool_c
      procedure :: get_eri_t1                => get_eri_t1_t1_eri_tool_c
      procedure :: get_g_pqrs_t1             => get_g_pqrs_t1_t1_eri_tool_c
      procedure :: construct_g_pqrs_t1       => construct_g_pqrs_t1_t1_eri_tool_c
      procedure :: construct_g_t1_from_mo    => construct_g_t1_from_mo_t1_eri_tool
!
!     Get and construct routines for packed g_pqrs in T1 basis
!
      procedure :: get_eri_t1_packed_mem      => get_eri_t1_packed_mem_t1_eri_tool_c
      procedure :: get_eri_t1_packed          => get_eri_t1_packed_t1_eri_tool_c
      procedure :: get_g_pqrs_t1_packed       => get_g_pqrs_t1_packed_t1_eri_tool_c
      procedure :: construct_g_pqrs_t1_packed => construct_g_pqrs_t1_packed_t1_eri_tool_c
!
!     Construction of C1 Cholesky vectors from T1 Cholesky
!
      procedure :: get_eri_c1_mem            => get_eri_c1_mem_t1_eri_tool_c
      procedure :: get_eri_c1                => get_eri_c1_t1_eri_tool_c
      procedure :: construct_g_pqrs_c1       => construct_g_pqrs_c1_t1_eri_tool_c
!
      procedure :: construct_cholesky_c1_mem => construct_cholesky_c1_mem_t1_eri_tool_c
      procedure :: construct_cholesky_c1     => construct_cholesky_c1_t1_eri_tool_c
      procedure :: construct_cholesky_oo_c1  => construct_cholesky_oo_c1_t1_eri_tool_c
      procedure :: construct_cholesky_vo_c1  => construct_cholesky_vo_c1_t1_eri_tool_c
      procedure :: construct_cholesky_vv_c1  => construct_cholesky_vv_c1_t1_eri_tool_c
!
!     Various routines
!
      procedure :: set_pointer_t1            => set_pointer_t1_t1_eri_tool_c
      procedure :: L_pointer_setup_t1        => L_pointer_setup_t1_t1_eri_tool_c
!
   end type t1_eri_tool_c
!
!
   interface t1_eri_tool_c
!
      procedure :: new_t1_eri_tool_c
      procedure :: new_t1_eri_tool_from_mo_c
      procedure :: new_t1_eri_tool_from_t1_c
      procedure :: new_t1_eri_tool_from_mo_r
      procedure :: new_t1_eri_tool_from_t1_r
!
   end interface t1_eri_tool_c
!
!
contains
!
!
   function new_t1_eri_tool_c(n_o, n_v, n_J, need_g_abcd) result(eri)
!!
!!    New T1 integral tool
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Creates the integral tool.
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
      type(t1_eri_tool_c) :: eri
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
      eri%cholesky_t1 = direct_stream_file('cholesky_t1', eri%n_J)
!
!     Default: Cholesky in memory, integral matrix in memory
!              depends on amount of memory and "need_g_abcd"
!
      if (eri%room_for_cholesky(4*dp)) then
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
        eri%t1_eri_mem = .true.
!
      else
!
        eri%t1_eri_mem = .false.
!
      endif
!
      eri%mo_eri_mem = .false.
!
!     Read non-default settings
!
      call eri%read_settings()
!
!     Print settings
!
      call eri%print_settings()
!
   end function new_t1_eri_tool_c
!
!
   function new_t1_eri_tool_from_mo_c(template) result(eri)
!!
!!    New T1 integral tool from MO complex
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Creates the integral tool based on template.
!!
      implicit none
!
      type(t1_eri_tool_c)  :: eri
!
      type(mo_eri_tool_c) :: template
!
      eri%n_J  = template%n_J
      eri%n_o  = template%n_o
      eri%n_v  = template%n_v
      eri%n_mo = template%n_mo
!
      eri%cholesky_mo = direct_stream_file('cholesky_mo', eri%n_J)
      eri%cholesky_t1 = direct_stream_file('cholesky_t1', eri%n_J)
!
!     Memory logicals - Cholesky in memory (or file)? ERI in memory (or file)?
!
      eri%cholesky_mem = template%cholesky_mem
      eri%mo_eri_mem   = template%mo_eri_mem
      eri%t1_eri_mem   = template%mo_eri_mem
!
!     Print settings
!
      call eri%print_settings()
!
   end function new_t1_eri_tool_from_mo_c
!
!
   function new_t1_eri_tool_from_t1_c(template) result(eri)
!!
!!    New T1 integral tool from T1 complex
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Creates the integral tool based on template.
!!
      implicit none
!
      type(t1_eri_tool_c)  :: eri
!
      type(t1_eri_tool_c) :: template
!
      eri%n_J  = template%n_J
      eri%n_o  = template%n_o
      eri%n_v  = template%n_v
      eri%n_mo = template%n_mo
!
      eri%cholesky_mo = direct_stream_file('cholesky_mo', eri%n_J)
      eri%cholesky_t1 = direct_stream_file('cholesky_t1', eri%n_J)
!
!     Memory logicals - Cholesky in memory (or file)? ERI in memory (or file)?
!
      eri%cholesky_mem = template%cholesky_mem
      eri%mo_eri_mem   = template%mo_eri_mem
      eri%t1_eri_mem   = template%t1_eri_mem
!
!     Print settings
!
      call eri%print_settings()
!
   end function new_t1_eri_tool_from_t1_c
!
!
   function new_t1_eri_tool_from_mo_r(template) result(eri)
!!
!!    New T1 integral tool from MO real
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Creates the integral tool based on template.
!!
      implicit none
!
      type(t1_eri_tool_c)  :: eri
!
      type(mo_eri_tool)    :: template
!
      eri%n_J  = template%n_J
      eri%n_o  = template%n_o
      eri%n_v  = template%n_v
      eri%n_mo = template%n_mo
!
      eri%cholesky_mo = direct_stream_file('cholesky_mo', eri%n_J)
      eri%cholesky_t1 = direct_stream_file('cholesky_t1', eri%n_J)
!
!     Memory logicals - Cholesky in memory (or file)? ERI in memory (or file)?
!
      eri%cholesky_mem = template%cholesky_mem
      eri%mo_eri_mem   = template%mo_eri_mem
      eri%t1_eri_mem   = template%mo_eri_mem
!
!     Print settings
!
      call eri%print_settings()
!
   end function new_t1_eri_tool_from_mo_r
!
!
   function new_t1_eri_tool_from_t1_r(template) result(eri)
!!
!!    New T1 integral tool from T1 real
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Creates the integral tool based on template.
!!
      implicit none
!
      type(t1_eri_tool_c)  :: eri
!
      type(t1_eri_tool)    :: template
!
      eri%n_J  = template%n_J
      eri%n_o  = template%n_o
      eri%n_v  = template%n_v
      eri%n_mo = template%n_mo
!
      eri%cholesky_mo = direct_stream_file('cholesky_mo', eri%n_J)
      eri%cholesky_t1 = direct_stream_file('cholesky_t1', eri%n_J)
!
!     Memory logicals - Cholesky in memory (or file)? ERI in memory (or file)?
!
      eri%cholesky_mem = template%cholesky_mem
      eri%mo_eri_mem   = template%mo_eri_mem
      eri%t1_eri_mem   = template%t1_eri_mem
!
!     Print settings
!
      call eri%print_settings()
!
   end function new_t1_eri_tool_from_t1_r
!
!
   subroutine print_settings_t1_eri_tool_c(eri)
!!
!!    Print settings
!!    Written by Eirik F. Kjønstad, Dec 2019
!!
      implicit none
!
      class(t1_eri_tool_c) :: eri
!
      call eri%mo_eri_tool_c%print_settings()
!
      call output%printf('m', &
                         'T1 ERI matrix in memory:    (l0)', &
                         logs=[eri%t1_eri_mem],              &
                         fs='(t6,a)')
!
   end subroutine print_settings_t1_eri_tool_c
!
!
   subroutine read_settings_t1_eri_tool_c(eri)
!!
!!    Read settings
!!    Written by Eirik F. Kjønstad, Dec 2019
!!
      implicit none
!
      class(t1_eri_tool_c) :: eri
!
      call eri%read_general_settings(eri%t1_eri_mem)
!
      if (input%is_keyword_present('t1 eri in memory', 'integrals')) then
!
            eri%t1_eri_mem = .true.
!
      endif
!
   end subroutine read_settings_t1_eri_tool_c
!
!
   subroutine initialize_t1_eri_tool_c(eri)
!!
!!    Initialize
!!    Written by Rolf H. Myhre, Aug 2020
!!
!!    Allocates arrays for Cholesky vectors and the ERIs
!!    if they are to be stored in memory.
!!
      implicit none
!
      class(t1_eri_tool_c) :: eri
!
      call eri%mo_eri_tool_c%initialize()
!
!     Allocate Cholesky arrays
!
      if (eri%cholesky_mem) then
!
         call mem%alloc(eri%L_J_oo_t1, eri%n_J, eri%n_o, eri%n_o)
         call mem%alloc(eri%L_J_ov_t1, eri%n_J, eri%n_o, eri%n_v)
         call mem%alloc(eri%L_J_vo_t1, eri%n_J, eri%n_v, eri%n_o)
         call mem%alloc(eri%L_J_vv_t1, eri%n_J, eri%n_v, eri%n_v)
!
      endif
!
!     Allocate integrals
!
      if (eri%t1_eri_mem) then
         call mem%alloc(eri%g_pqrs_t1, eri%n_mo, eri%n_mo, &
                                       eri%n_mo, eri%n_mo)
      endif
!
   end subroutine initialize_t1_eri_tool_c
!
!
   subroutine copy_from_mo_t1_eri_tool_c(eri, template)
!!
!!    Copy from MO
!!    Written by Rolf H. Myhre, Aug 2020
!!
!!    Copies MO Cholesky and ERIs into MO and T1 Cholesky and ERIs
!!    from complex MO template if in memory
!!
      implicit none
!
      class(t1_eri_tool_c) :: eri
!
      class(mo_eri_tool_c) :: template
!
!     Copy Cholesky vectors from template - if present
!
      call eri%mo_eri_tool_c%copy_from_mo(template)
!
      if (eri%cholesky_mem .and. template%cholesky_mem) then
!
         call zcopy(eri%n_J*eri%n_o*eri%n_o, &
                    template%L_J_oo_mo,  &
                    1, eri%L_J_oo_t1, 1)
!
         call zcopy(eri%n_J*eri%n_o*eri%n_v, &
                    template%L_J_ov_mo,  &
                    1, eri%L_J_ov_t1, 1)
!
         call zcopy(eri%n_J*eri%n_v*eri%n_o, &
                    template%L_J_vo_mo,  &
                    1, eri%L_J_vo_t1, 1)
!
         call zcopy(eri%n_J*eri%n_v*eri%n_v, &
                    template%L_J_vv_mo,  &
                    1, eri%L_J_vv_t1, 1)
!
      endif
!
!
!     Copy integrals from template
!
      if (eri%t1_eri_mem .and. template%mo_eri_mem) then
         call zcopy(eri%n_mo**4,        &
                    template%g_pqrs_mo, &
                    1,                  &
                    eri%g_pqrs_t1,      &
                    1)
      endif
!
   end subroutine copy_from_mo_t1_eri_tool_c
!
!
   subroutine copy_from_real_mo_t1_eri_tool_c(eri, template)
!!
!!    Copy from real MO
!!    Written by Rolf H. Myhre, Aug 2020
!!
!!    Copies MO Cholesky and ERIs into MO and T1 Cholesky and ERIs
!!    from real MO template if in memory
!!
      implicit none
!
      class(t1_eri_tool_c) :: eri
!
      class(mo_eri_tool)   :: template
!
!     Copy Cholesky vectors from template - if present
!
      call eri%mo_eri_tool_c%copy_from_real_mo(template)
!
      if (eri%cholesky_mem .and. template%cholesky_mem) then
!
         eri%L_J_oo_t1 = cmplx(template%L_J_oo_mo, zero, dp)
         eri%L_J_vo_t1 = cmplx(template%L_J_vo_mo, zero, dp)
         eri%L_J_ov_t1 = cmplx(template%L_J_ov_mo, zero, dp)
         eri%L_J_vv_t1 = cmplx(template%L_J_vv_mo, zero, dp)
!
      endif
!
!
!     Copy integrals from template - if present
!
      if (eri%t1_eri_mem .and. template%mo_eri_mem) then
         eri%g_pqrs_t1 = cmplx(template%g_pqrs_mo, zero, dp)
      endif
!
   end subroutine copy_from_real_mo_t1_eri_tool_c
!
!
   subroutine copy_from_t1_t1_eri_tool_c(eri, template)
!!
!!    Copy from T1
!!    Written by Rolf H. Myhre, Aug 2020
!!
!!    Copies Cholesky vectors and ERIs from complex T1 template if in memory
!!
      implicit none
!
      class(t1_eri_tool_c) :: eri
!
      class(t1_eri_tool_c) :: template
!
!     Copy Cholesky vectors from template - if present
!
      call eri%mo_eri_tool_c%copy_from_mo(template)
!
      if (eri%cholesky_mem .and. template%cholesky_mem) then
!
         call zcopy(eri%n_J*eri%n_o*eri%n_o, &
                    template%L_J_oo_t1,  &
                    1, eri%L_J_oo_t1, 1)
!
         call zcopy(eri%n_J*eri%n_o*eri%n_v, &
                    template%L_J_ov_t1,  &
                    1, eri%L_J_ov_t1, 1)
!
         call zcopy(eri%n_J*eri%n_v*eri%n_o, &
                    template%L_J_vo_t1,  &
                    1, eri%L_J_vo_t1, 1)
!
         call zcopy(eri%n_J*eri%n_v*eri%n_v, &
                    template%L_J_vv_t1,  &
                    1, eri%L_J_vv_t1, 1)
!
      endif
!
!
!     Copy integrals from template - if present
!
      if (eri%t1_eri_mem .and. template%t1_eri_mem) then
         call zcopy(eri%n_mo**4,        &
                    template%g_pqrs_t1, &
                    1,                  &
                    eri%g_pqrs_t1,      &
                    1)
      endif
!
   end subroutine copy_from_t1_t1_eri_tool_c
!
!
   subroutine copy_from_real_t1_t1_eri_tool_c(eri, template)
!!
!!    Copy from real T1
!!    Written by Rolf H. Myhre, Aug 2020
!!
!!    Copies Cholesky vectors and ERIs from real T1 template if in memory
!!
      implicit none
!
      class(t1_eri_tool_c) :: eri
!
      class(t1_eri_tool)   :: template
!
!     Copy Cholesky vectors from template - if present
!
      call eri%mo_eri_tool_c%copy_from_real_mo(template)
!
      if (eri%cholesky_mem .and. template%cholesky_mem) then
!
         eri%L_J_oo_t1 = cmplx(template%L_J_oo_t1, zero, dp)
         eri%L_J_vo_t1 = cmplx(template%L_J_vo_t1, zero, dp)
         eri%L_J_ov_t1 = cmplx(template%L_J_ov_t1, zero, dp)
         eri%L_J_vv_t1 = cmplx(template%L_J_vv_t1, zero, dp)
!
      endif
!
!
!     Copy integrals from template - if present
!
      if (eri%t1_eri_mem .and. template%t1_eri_mem) then
         eri%g_pqrs_t1 = cmplx(template%g_pqrs_t1, zero, dp)
      endif
!
   end subroutine copy_from_real_t1_t1_eri_tool_c
!
!
   subroutine set_cholesky_t1_t1_eri_tool_c(eri, L_J_pq, first_p, last_p, first_q, last_q)
!!
!!    Set Cholesky T1
!!    Written by Rolf H. Myhre, May 2020
!!
!!    based on routines by Sarai D. Folkestad and Eirik F. Kjønstad
!!
!!    Saves parts of the Cholesky T1 vectors either to file or in memory.
!!
!!    L_J_pq:           Cholesky vector L_pq^J ordered as J x pq
!!
!!    first_p, last_p:  MO index range for p in L_J_pq
!!    first_q, last_q:  MO index range for q in L_J_pq
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      complex(dp), dimension(eri%n_J, first_p:last_p, first_q:last_q), intent(in) :: L_J_pq
!
      if (last_p .lt. first_p) call output%error_msg('last_p .lt. first_p in set_cholesky_t1')
      if (last_q .lt. first_q) call output%error_msg('last_q .lt. first_q in set_cholesky_t1')
!
      if (eri%cholesky_mem) then
!
         call eri%copy_in_cholesky(L_J_pq, first_p, last_p, first_q, last_q, &
                                         eri%L_J_oo_t1, &
                                         eri%L_J_vo_t1, &
                                         eri%L_J_ov_t1, &
                                         eri%L_J_vv_t1)
      else
!
         call eri%cholesky_t1%open_('write')
!
         call eri%write_cholesky(L_J_pq, &
                                       first_p, last_p, &
                                       first_q, last_q, &
                                       eri%cholesky_t1)
!
         call eri%cholesky_t1%close_()
!
      endif
!
   end subroutine set_cholesky_t1_t1_eri_tool_c
!
!
   subroutine get_cholesky_t1_t1_eri_tool_c(eri, L_J_pq, first_p, last_p, first_q, last_q)
!!
!!    Get Cholesky T1
!!    Written by Rolf H. Myhre, May 2020
!!
!!    Saves parts of the Cholesky T1 vectors either to file or in memory.
!!
!!    L_J_pq:           Cholesky vector L_pq^J ordered as J x pq
!!
!!    first_p, last_p:  MO index range for p in L_J_pq
!!    first_q, last_q:  MO index range for q in L_J_pq
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      complex(dp), dimension(eri%n_J, first_p:last_p, first_q:last_q), intent(out) :: L_J_pq
!
      if (last_p .lt. first_p) call output%error_msg('last_p .lt. first_p in get_cholesky_t1')
      if (last_q .lt. first_q) call output%error_msg('last_q .lt. first_q in get_cholesky_t1')
!
      if (eri%cholesky_mem) then
!
         call eri%copy_out_cholesky(L_J_pq, first_p, last_p, first_q, last_q, &
                                          eri%L_J_oo_t1, &
                                          eri%L_J_vo_t1, &
                                          eri%L_J_ov_t1, &
                                          eri%L_J_vv_t1)
      else
!
         call eri%cholesky_t1%open_('read')
!
         call eri%read_cholesky(L_J_pq, &
                                      first_p, last_p, &
                                      first_q, last_q, &
                                      eri%cholesky_t1)
!
         call eri%cholesky_t1%close_()
!
      endif
!
   end subroutine get_cholesky_t1_t1_eri_tool_c
!
!
   subroutine set_t1_to_mo_t1_eri_tool_c(eri)
!!
!!    Set T1 to MO
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    Sets the T1 Cholesky vector equal to the MO vectors
!!    and constructs the MO ERIs if in memory
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      type(batching_index) :: batcher
      integer :: batch
!
      complex(dp), dimension(:), allocatable :: array
!
      type(timings) :: copy_cholesky_timer, construct_eri_timer
!
      copy_cholesky_timer = timings('Copy Cholesky MO to T1', pl='normal')
      call copy_cholesky_timer%turn_on()
!
      if (eri%cholesky_mem) then !Copy everything
         call zcopy(eri%n_J*eri%n_o*eri%n_o, &
                    eri%L_J_oo_mo, 1, &
                    eri%L_J_oo_t1, 1)
!
         call zcopy(eri%n_J*eri%n_v*eri%n_o, &
                    eri%L_J_vo_mo, 1, &
                    eri%L_J_vo_t1, 1)
!
         call zcopy(eri%n_J*eri%n_o*eri%n_v, &
                    eri%L_J_ov_mo, 1, &
                    eri%L_J_ov_t1, 1)
!
         call zcopy(eri%n_J*eri%n_v*eri%n_v, &
                    eri%L_J_vv_mo, 1, &
                    eri%L_J_vv_t1, 1)
!
      else !Read and write
!
         batcher = batching_index(eri%n_mo**2)
         call mem%batch_setup(batcher, 0, eri%n_J)
         call mem%alloc(array, batcher%max_length*eri%n_J)
!
         call eri%cholesky_mo%open_('read')
         call eri%cholesky_t1%open_('write')
!
         do batch = 1, batcher%num_batches
!
            call batcher%determine_limits(batch)
            call eri%cholesky_mo%read_(array, batcher%first, batcher%last)
            call eri%cholesky_t1%write_(array, batcher%first, batcher%last)
!
         enddo
!
         call mem%dealloc(array, batcher%max_length*eri%n_J)
         call eri%cholesky_mo%close_()
         call eri%cholesky_t1%close_()
!
      endif
      call copy_cholesky_timer%turn_off()
!
      if (eri%t1_eri_mem) then
!
         construct_eri_timer = timings('Construct eri integrals', pl='normal')
         call construct_eri_timer%turn_on()
!
         call eri%construct_g_pqrs_t1(eri%g_pqrs_t1, &
                                      1, eri%n_mo,   &
                                      1, eri%n_mo,   &
                                      1, eri%n_mo,   &
                                      1, eri%n_mo,   &
                                      one_complex, zero_complex)
!
         call construct_eri_timer%turn_off()
!
      endif
!
   end subroutine set_t1_to_mo_t1_eri_tool_c
!
!
   subroutine update_t1_integrals_t1_eri_tool_c(eri, t1)
!!
!!    Update T1 integrals
!!    Written by Eirik F. Kjønstad, Dec 2019
!!
!!    Updates the T1 integrals by:
!!
!!       1. Constructing and saving (in memory or file) the Cholesky
!!          T1 transformed vectors
!!
!!       2. If ERI in memory, constructs and saves T1 g_pqrs
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      complex(dp), dimension(eri%n_v, eri%n_o), intent(in) :: t1
!
      type(timings) :: update_t1_cholesky_timer, construct_eri_timer
!
      update_t1_cholesky_timer = timings('Update t1 cholesky', pl='normal')
      call update_t1_cholesky_timer%turn_on()
!
      call eri%construct_t1_cholesky(t1)
!
      call update_t1_cholesky_timer%turn_off()
!
      if (eri%t1_eri_mem) then
!
         construct_eri_timer = timings('Construct eri integrals', pl='normal')
         call construct_eri_timer%turn_on()
!
         if (eri%mo_eri_mem) then
!
            call eri%construct_g_t1_from_mo(t1)
!
         else
!
            call eri%construct_g_pqrs_t1(eri%g_pqrs_t1, &
                                         1, eri%n_mo,   &
                                         1, eri%n_mo,   &
                                         1, eri%n_mo,   &
                                         1, eri%n_mo,   &
                                         one_complex, zero_complex)
!
         endif
!
         call construct_eri_timer%turn_off()
!
      endif
!
   end subroutine update_t1_integrals_t1_eri_tool_c
!
!
   subroutine construct_t1_cholesky_t1_eri_tool_c(eri, t1)
!!
!!    Construct T1 Cholesky
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    based on routines by Sarai D. Folkestad and Eirik F. Kjønstad
!!
!!    Wrapper for the constructors of the different Cholesky T1 blocks
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      complex(dp), dimension(eri%n_v, eri%n_o), intent(in) ::t1
!
!     occupied-occupied block
!
      call eri%construct_cholesky_t1_oo(t1)
!
!     virtual-occupied block
!
      call eri%construct_cholesky_t1_vo(t1)
!
!     virtual-virtual block
!
      call eri%construct_cholesky_t1_vv(t1)
!
   end subroutine construct_t1_cholesky_t1_eri_tool_c
!
!
   subroutine construct_cholesky_t1_oo_t1_eri_tool_c(eri, t1)
!!
!!    Construct Cholesky T1 oo
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    based on routines by Sarai D. Folkestad and Eirik F. Kjønstad
!!
!!    Computes
!!
!!       L_J_ij_T1 = L_J_ij_MO + sum_a L_J_ia_MO*t_aj,
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      complex(dp), dimension(eri%n_v, eri%n_o), intent(in) :: t1
!
      complex(dp), dimension(:), allocatable :: L_J_oo
      complex(dp), dimension(:), allocatable :: L_J_ox
!
      type(batching_index) :: batch_o
!
      integer :: o_batch
!
      if (eri%cholesky_mem) then
!
!        L_J_ij_t1 = L_J_ij_mo
         call zcopy(eri%n_J*eri%n_o**2, eri%L_J_oo_mo, 1, eri%L_J_oo_t1, 1)
!
!        L_J_ij_t1 += sum_b L_J_ib_mo t_bj
         call zgemm('N', 'N',                          &
                    eri%n_J*eri%n_o, eri%n_o, eri%n_v, &
                    one_complex,                       &
                    eri%L_J_ov_mo, eri%n_J*eri%n_o,    &
                    t1, eri%n_v,                       &
                    one_complex,                       &
                    eri%L_J_oo_t1, eri%n_J*eri%n_o)
!
      else
!
         batch_o = batching_index(eri%n_o)
         call mem%batch_setup(batch_o, 0, eri%n_J*(eri%n_o*max(eri%n_v, eri%n_o) + eri%n_o**2))
!
         call mem%alloc(L_J_oo, eri%n_J*batch_o%max_length*eri%n_o)
         call mem%alloc(L_J_ox, eri%n_J*batch_o%max_length*max(eri%n_v, eri%n_o))
!
         do o_batch = 1,batch_o%num_batches
!
            call batch_o%determine_limits(o_batch)
!
!           L_J_ij_t1 = L_J_ij_mo
            call eri%get_cholesky_mo(L_J_oo, batch_o%first, batch_o%last, 1, eri%n_o)
!
!           L_J_ij_t1 += sum_b L_J_ib_mo t_bj
            call eri%get_cholesky_mo(L_J_ox, batch_o%first, batch_o%last, eri%n_o+1, eri%n_mo)
!
            call zgemm('N', 'N',                                 &
                       eri%n_J*batch_o%length, eri%n_o, eri%n_v, &
                       one_complex,                              &
                       L_J_ox, eri%n_J*batch_o%length,           &
                       t1, eri%n_v,                              &
                       one_complex,                              &
                       L_J_oo, eri%n_J*batch_o%length)
!
            call eri%set_cholesky_t1(L_J_oo, batch_o%first, batch_o%last, 1, eri%n_o)
!
         enddo
!
         call mem%dealloc(L_J_oo, eri%n_J*batch_o%max_length*eri%n_o)
         call mem%dealloc(L_J_ox, eri%n_J*batch_o%max_length*max(eri%n_v, eri%n_o))
!
      end if
!
   end subroutine construct_cholesky_t1_oo_t1_eri_tool_c
!
!
   subroutine construct_cholesky_t1_vo_t1_eri_tool_c(eri, t1)
!!
!!    Construct Cholesky T1 vo
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    based on routines by Sarai D. Folkestad and Eirik F. Kjønstad
!!
!!    Computes
!!
!!       L_J_ai_T1 = L_J_ai_MO + sum_b L_J_ab_MO*t_bi
!!                             - sum_j t_aj*(L_J_ji_MO + sum_b L_J_jb_MO*t_bi)
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      complex(dp), dimension(eri%n_v, eri%n_o), intent(in) :: t1
!
      complex(dp), dimension(:), allocatable, target :: L_J_ij
      complex(dp), dimension(:), allocatable, target :: L_J_ji
      complex(dp), dimension(:), allocatable, target :: L_J_ai
      complex(dp), dimension(:), allocatable, target :: L_J_ab
!
      complex(dp), dimension(:,:,:), allocatable :: L_J_ia
!
      complex(dp), dimension(:,:,:), pointer :: L_J_ij_p
      complex(dp), dimension(:,:,:), pointer :: L_J_ji_p
      complex(dp), dimension(:,:,:), pointer :: L_J_ia_p
      complex(dp), dimension(:,:,:), pointer :: L_J_ai_p
      complex(dp), dimension(:,:,:), pointer :: L_J_ab_p
!
      type(batching_index) :: batch_v, batch_o
!
      integer :: v_batch, o_batch, i, a, J, req
!
      if (eri%cholesky_mem) then
!
!        L_J_ai_t1 = L_J_ai_mo
         call zcopy(eri%n_J*eri%n_v*eri%n_o, eri%L_J_vo_mo, 1, eri%L_J_vo_t1, 1)
!
!        L_J_ai_t1 += sum_b L_J_ab_mo t_bi
         call zgemm('N', 'N',                          &
                    eri%n_J*eri%n_v, eri%n_o, eri%n_v, &
                    one_complex,                       &
                    eri%L_J_vv_mo, eri%n_J*eri%n_v,    &
                    t1, eri%n_v,                       &
                    one_complex,                       &
                    eri%L_J_vo_t1, eri%n_J*eri%n_v)
!
         batch_o = batching_index(eri%n_o)
         batch_v = batching_index(eri%n_v)
         call mem%batch_setup(batch_o, batch_v, 0, 2*eri%n_J*eri%n_o, eri%n_J*eri%n_o, 0)
!
         call mem%alloc(L_J_ij, eri%n_J*eri%n_o*batch_o%max_length)
         call mem%alloc(L_J_ji, eri%n_J*batch_o%max_length*eri%n_o)
         call mem%alloc(L_J_ia, eri%n_J, eri%n_o, batch_v%max_length)
!
         do o_batch = 1, batch_o%num_batches
!
            call batch_o%determine_limits(o_batch)
!
            L_J_ji_p(1:eri%n_J, 1:batch_o%length, 1:eri%n_o) &
                 => L_J_ji(1:eri%n_J*batch_o%length*eri%n_o)
!
            L_J_ij_p(1:eri%n_J, 1:eri%n_o, 1:batch_o%length) &
                 => L_J_ij(1:eri%n_J*batch_o%length*eri%n_o)
!
!           L'_J_ji = sum_b L_J_jb_mo t_bi
            call zgemm('N', 'N', &
                       eri%n_J*batch_o%length, eri%n_o, eri%n_v,          &
                       one_complex,                                       &
                       eri%L_J_ov_mo(:,batch_o%first,1), eri%n_J*eri%n_o, &
                       t1, eri%n_v,                                       &
                       zero_complex,                                      &
                       L_J_ji, eri%n_J*batch_o%length)
!
            call sort_123_to_132(L_J_ji_p, L_J_ij_p, eri%n_J, batch_o%length, eri%n_o)
!
!           L'_J_ji += L_J_ji_mo
            call zaxpy(eri%n_J*eri%n_o*batch_o%length, one_complex, &
                       eri%L_J_oo_mo(:,:,batch_o%first), 1, L_J_ij,1)
!
            do v_batch = 1,batch_v%num_batches
!
               call batch_v%determine_limits(v_batch)
!
!              L_J_ai_t1 -= sum_j t_aj L'_J_ji_mo
               call zgemm('N', 'T',                                        &
                          eri%n_J*eri%n_o, batch_v%length, batch_o%length, &
                          one_complex,                                     &
                          L_J_ij, eri%n_J*eri%n_o,                         &
                          t1(batch_v%first:, batch_o%first), eri%n_v,      &
                          zero_complex,                                    &
                          L_J_ia, eri%n_J*eri%n_o)
!
!$omp parallel do private(i, a, J)
               do i = 1, eri%n_o
                  do a = 1, batch_v%length
                     do J = 1, eri%n_J
                        eri%L_J_vo_t1(J, batch_v%first+a-1, i) = &
                        eri%L_J_vo_t1(J, batch_v%first+a-1, i) - L_J_ia(J,i,a)
                     enddo
                  enddo
               enddo
!$omp end parallel do
!
            enddo
         enddo
!
         call mem%dealloc(L_J_ia, eri%n_J, eri%n_o, batch_v%max_length)
         call mem%dealloc(L_J_ji, eri%n_J*batch_o%max_length*eri%n_o)
         call mem%dealloc(L_J_ij, eri%n_J*eri%n_o*batch_o%max_length)
!
      else
!
         batch_o = batching_index(eri%n_o)
         batch_v = batching_index(eri%n_v)
         req = eri%n_J*(eri%n_o + max(eri%n_v,eri%n_o))
         call mem%batch_setup(batch_o, batch_v, 0, req, req, 0)
!
         call mem%alloc(L_J_ij, eri%n_J*batch_o%max_length*max(eri%n_v, eri%n_o))
         call mem%alloc(L_J_ji, eri%n_J*batch_o%max_length*eri%n_o)
!
         call mem%alloc(L_J_ai, eri%n_J*batch_v%max_length*eri%n_o)
         call mem%alloc(L_J_ab, eri%n_J*batch_v%max_length*max(eri%n_v, eri%n_o))
!
         do v_batch = 1, batch_v%num_batches
!
            call batch_v%determine_limits(v_batch)
!
            L_J_ai_p(1:eri%n_J, 1:batch_v%length, 1:eri%n_o) => L_J_ai
            L_J_ab_p(1:eri%n_J, 1:batch_v%length, 1:eri%n_v) => L_J_ab
!
!           L_J_ai_t1 = L_J_ai_mo
            call eri%get_cholesky_mo(L_J_ai_p, &
                                     eri%n_o+batch_v%first, eri%n_o+batch_v%last, &
                                     1, eri%n_o)
!
!           L_J_ai_t1 += sum_b L_J_ab_mo t_bi
            call eri%get_cholesky_mo(L_J_ab_p, &
                                     eri%n_o+batch_v%first, eri%n_o+batch_v%last, &
                                     eri%n_o+1, eri%n_mo)
!
            call zgemm('N', 'N',                                 &
                       eri%n_J*batch_v%length, eri%n_o, eri%n_v, &
                       one_complex,                              &
                       L_J_ab_p, eri%n_J*batch_v%length,         &
                       t1, eri%n_v,                              &
                       one_complex,                              &
                       L_J_ai_p, eri%n_J*batch_v%length)
!
            do o_batch = 1, batch_o%num_batches
!
               call batch_o%determine_limits(o_batch)
!
               L_J_ji_p(1:eri%n_J, 1:batch_o%length, 1:eri%n_o) => L_J_ji
               L_J_ia_p(1:eri%n_J, 1:batch_o%length, 1:eri%n_v) => L_J_ij
!
!              L'_J_ji = sum_b L_J_jb_mo t_bi
               call eri%get_cholesky_mo(L_J_ia_p, batch_o%first, batch_o%last, eri%n_o+1, eri%n_mo)
!
               call zgemm('N', 'N',                                 &
                          eri%n_J*batch_o%length, eri%n_o, eri%n_v, &
                          one_complex,                              &
                          L_J_ia_p, eri%n_J*batch_o%length,         &
                          t1, eri%n_v,                              &
                          zero_complex,                             &
                          L_J_ji_p, eri%n_J*batch_o%length)
!
!              L'_J_ji += L_J_ji_mo
               L_J_ij_p(1:eri%n_J, 1:eri%n_o, 1:batch_o%length) => L_J_ij
               call eri%get_cholesky_mo(L_J_ij_p, 1, eri%n_o, batch_o%first, batch_o%last)
!
               call add_132_to_123(one_complex, L_J_ji_p, L_J_ij_p, &
                                   eri%n_J, eri%n_o, batch_o%length)
!
!              L_J_ai_t1 -= sum_j t_aj L'_J_ji
               L_J_ia_p(1:eri%n_J, 1:eri%n_o, 1:batch_v%length) => L_J_ab
!
               call zgemm('N', 'T',                                        &
                          eri%n_J*eri%n_o, batch_v%length, batch_o%length, &
                          -one_complex,                                    &
                          L_J_ij_p, eri%n_J*eri%n_o,                       &
                          t1(batch_v%first:, batch_o%first), eri%n_v,      &
                          zero_complex,                                    &
                          L_J_ia_p, eri%n_J*eri%n_o)
!
               call add_132_to_123(one_complex, L_J_ia_p, L_J_ai_p, &
                                   eri%n_J, batch_v%length, eri%n_o)
!
            enddo
!
            call eri%set_cholesky_t1(L_J_ai_p, &
                                     eri%n_o+batch_v%first, eri%n_o+batch_v%last, 1, eri%n_o)
!
         enddo
!
         call mem%dealloc(L_J_ab, eri%n_J*batch_v%max_length*max(eri%n_v, eri%n_o))
         call mem%dealloc(L_J_ai, eri%n_J*batch_v%max_length*eri%n_o)
!
         call mem%dealloc(L_J_ji, eri%n_J*batch_o%max_length*eri%n_o)
         call mem%dealloc(L_J_ij, eri%n_J*batch_o%max_length*max(eri%n_v, eri%n_o))
!
      endif
!
   end subroutine construct_cholesky_t1_vo_t1_eri_tool_c
!
!
   subroutine construct_cholesky_t1_vv_t1_eri_tool_c(eri, t1)
!!
!!    Construct Cholesky T1 vv
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    based on routines by Sarai D. Folkestad and Eirik F. Kjønstad
!!
!!    Computes
!!
!!       L_J_ab_T1 = L_J_ab_MO - sum_j t_aj*L_J_jb_MO
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      complex(dp), dimension(eri%n_v, eri%n_o), intent(in) :: t1
!
      complex(dp), dimension(:,:,:), allocatable :: L_J_vv
      complex(dp), dimension(:,:,:), allocatable, target :: L_J_xv
!
      type(batching_index) :: batch_v
!
      integer :: v_batch
!
      if (eri%cholesky_mem) then
!
!        L_J_ab_t1 = L_J_ab_mo
         call zcopy(eri%n_J*eri%n_v**2, eri%L_J_vv_mo, 1, eri%L_J_vv_t1, 1)
!
         batch_v = batching_index(eri%n_v)
         call mem%batch_setup(batch_v, 0, eri%n_J*eri%n_v)
!
         call mem%alloc(L_J_vv, eri%n_J, eri%n_v, batch_v%max_length)
!
         do v_batch = 1, batch_v%num_batches
!
            call batch_v%determine_limits(v_batch)
!
!           L_J_ab_t1 -= sum_j t_aj L_J_jb_mo
            call zgemm('N', 'T',                                           &
                       eri%n_J*batch_v%length, eri%n_v, eri%n_o,           &
                       -one_complex,                                       &
                       eri%L_J_vo_mo(:,batch_v%first, 1), eri%n_J*eri%n_v, &
                       t1, eri%n_v,                                        &
                       zero_complex,                                       &
                       L_J_vv, eri%n_J*batch_v%length)
!
            call add_132_to_123(one_complex, L_J_vv, &
                                eri%L_J_vv_t1(:,:,batch_v%first:batch_v%last), &
                                eri%n_J, eri%n_v, batch_v%length)
!
         enddo
!
         call mem%dealloc(L_J_vv, eri%n_J, eri%n_v, batch_v%max_length)
!
      else
!
         batch_v = batching_index(eri%n_v)
         call mem%batch_setup(batch_v, 0, 2*eri%n_J*max(eri%n_v, eri%n_o))
!
         call mem%alloc(L_J_xv, eri%n_J, batch_v%max_length, max(eri%n_v, eri%n_o))
         call mem%alloc(L_J_vv, eri%n_J, batch_v%max_length, max(eri%n_v, eri%n_o))
!
         do v_batch = 1, batch_v%num_batches
!
            call batch_v%determine_limits(v_batch)
!
!           L_J_ab_t1 -= sum_j t_aj L_J_jb_mo
            call eri%get_cholesky_mo(L_J_xv, &
                                     eri%n_o + batch_v%first, eri%n_o + batch_v%last, 1, eri%n_o)
!
            call zgemm('N', 'T',                                 &
                       eri%n_J*batch_v%length, eri%n_v, eri%n_o, &
                       -one_complex,                             &
                       L_J_xv, eri%n_J*batch_v%length,           &
                       t1, eri%n_v,                              &
                       zero_complex,                             &
                       L_J_vv, eri%n_J*batch_v%length)
!
!           L_J_ab_t1 = L_J_ab_mo
            call eri%get_cholesky_mo(L_J_xv, eri%n_o + 1, eri%n_mo, &
                                     eri%n_o + batch_v%first, eri%n_o + batch_v%last)
!
            call add_132_to_123(one_complex, L_J_vv, L_J_xv, eri%n_J, eri%n_v, batch_v%length)
!
            call eri%set_cholesky_t1(L_J_xv, eri%n_o+1, eri%n_mo, &
                                     eri%n_o + batch_v%first, eri%n_o + batch_v%last)
!
         enddo
!
         call mem%dealloc(L_J_xv, eri%n_J, batch_v%max_length, max(eri%n_v, eri%n_o))
         call mem%dealloc(L_J_vv, eri%n_J, batch_v%max_length, max(eri%n_v, eri%n_o))
!
      end if
!
   end subroutine construct_cholesky_t1_vv_t1_eri_tool_c
!
!
   subroutine construct_g_t1_from_mo_t1_eri_tool(eri, t1)
!!
!!    Construct g T1 from MO
!!    Written by Rolf H. Myhre, Aug 2020
!!
!!    Construct g_pqrs in T1 basis from g_pqrs in MO basis
!!    In cases when both integrals are in memory,
!!    this is faster than constructing from L
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      complex(dp), dimension(eri%n_v, eri%n_o), intent(in) :: t1
!
      complex(dp), dimension(:,:,:,:), allocatable :: g_temp
!
      call mem%alloc(g_temp, eri%n_mo, eri%n_mo, eri%n_mo, eri%n_mo)
!
      call zcopy(eri%n_mo**4, eri%g_pqrs_mo, 1, eri%g_pqrs_t1, 1)
!
      call zgemm('N', 'N',                                    &
                 eri%n_mo**3, eri%n_o, eri%n_v,               &
                 one_complex,                                 &
                 eri%g_pqrs_mo(:,:,:,eri%n_o+1), eri%n_mo**3, &
                 t1, eri%n_v,                                 &
                 one_complex,                                 &
                 eri%g_pqrs_t1, eri%n_mo**3)
!
      call zcopy(eri%n_mo**4, eri%g_pqrs_t1, 1, g_temp, 1)
!
      call zgemm('N', 'N',                           &
                 eri%n_v, eri%n_mo**3, eri%n_o,      &
                 -one_complex,                       &
                 t1, eri%n_v,                        &
                 eri%g_pqrs_t1, eri%n_mo,            &
                 one_complex,                        &
                 g_temp(eri%n_o+1:,1,1,1), eri%n_mo)
!
      call sort_1234_to_3412(g_temp, eri%g_pqrs_t1, eri%n_mo, eri%n_mo, eri%n_mo, eri%n_mo)
!
      call zcopy(eri%n_mo**4, eri%g_pqrs_t1, 1, g_temp, 1)
!
      call zgemm('N', 'N',                                    &
                 eri%n_mo**3, eri%n_o, eri%n_v,               &
                 one_complex,                                 &
                 eri%g_pqrs_t1(:,:,:,eri%n_o+1), eri%n_mo**3, &
                 t1, eri%n_v,                                 &
                 one_complex,                                 &
                 g_temp, eri%n_mo**3)
!
      call zcopy(eri%n_mo**4, g_temp, 1, eri%g_pqrs_t1, 1)
!
      call zgemm('N', 'N',                           &
                 eri%n_v, eri%n_mo**3, eri%n_o,      &
                 -one_complex,                       &
                 t1, eri%n_v,                        &
                 g_temp, eri%n_mo,                   &
                 one_complex,                        &
                 eri%g_pqrs_t1(eri%n_o+1:,1,1,1), eri%n_mo)
!
      call mem%dealloc(g_temp, eri%n_mo, eri%n_mo, eri%n_mo, eri%n_mo)
!
   end subroutine construct_g_t1_from_mo_t1_eri_tool
!
!
   pure subroutine get_eri_t1_mem_t1_eri_tool_c(eri, string, req_pq, req_rs, &
                                                dim_p, dim_q, dim_r, dim_s,  &
                                                qp, sr)
!!
!!    Get ERI T1 mem
!!    Written by Rolf H. Myhre Jun 2020
!!
!!    Adds the memory usage for get ERI T1 depending on pq and rs to req_pq and req_rs
!!    Note that req_pq and req_rs must be initialized outside
!!    This routine is meant for estimates for the batching manager.
!!    If you are batching over index r, dim_r will typically be 1
!!    and the other dimensions should be full.
!!
!!    Temporary arrays needed for reordering, triggered by optional qp and sr
!!    are not allocated simultaneously in construct_eri, so this routine may
!!    overestimate total memory usage
!!
!!    See get_eri_t1 for additional documentation
!!
      implicit none
!
      class(t1_eri_tool_c), intent(in) :: eri
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
      if(.not. (eri%t1_eri_mem)) then
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
   end subroutine get_eri_t1_mem_t1_eri_tool_c
!
!
   subroutine get_eri_t1_t1_eri_tool_c(eri, string, g_pqrs, &
                                       first_p, last_p, first_q, last_q, &
                                       first_r, last_r, first_s, last_s, &
                                       alpha, beta, qp, sr, rspq)
!!
!!    Get ERI T1
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
      class(t1_eri_tool_c), intent(inout) :: eri
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
      call eri%get_g_pqrs_t1(g_pqrs, full_first_p, full_last_p, &
                                     full_first_q, full_last_q, &
                                     full_first_r, full_last_r, &
                                     full_first_s, full_last_s, &
                                     alpha_, beta_, qp, sr, rspq)
!
   end subroutine get_eri_t1_t1_eri_tool_c
!
!
   subroutine get_g_pqrs_t1_t1_eri_tool_c(eri, g_pqrs, first_p, last_p, first_q, last_q, &
                                                       first_r, last_r, first_s, last_s, &
                                                       alpha, beta, qp, sr, rspq)
!!
!!    Get g_pqrs T1
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
      class(t1_eri_tool_c), intent(inout) :: eri
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
      if (eri%t1_eri_mem) then
!
         call eri%copy_g_pqrs(g_pqrs, eri%g_pqrs_t1, &
                              first_p, last_p, first_q, last_q, &
                              first_r, last_r, first_s, last_s, &
                              alpha, beta, qp, sr, rspq)
!
      else
!
         call eri%construct_g_pqrs_t1(g_pqrs, first_p, last_p, first_q, last_q, &
                                              first_r, last_r, first_s, last_s, &
                                              alpha, beta, qp, sr, rspq)
!
      endif
!
   end subroutine get_g_pqrs_t1_t1_eri_tool_c
!
!
   subroutine construct_g_pqrs_t1_t1_eri_tool_c(eri, g_pqrs, &
                                           first_p, last_p, first_q, last_q, &
                                           first_r, last_r, first_s, last_s, &
                                           alpha, beta, qp, sr, rspq)
!!
!!    Construct g_pqrs T1
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
      class(t1_eri_tool_c), intent(inout) :: eri
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
      call eri%L_pointer_setup_t1(L_J_pq_p, switch_pq, &
                                  first_p, last_p, first_q, last_q, pq_alloced)
!
      if((first_p .eq. first_r) .and. (last_p .eq. last_r) .and. &
         (first_q .eq. first_s) .and. (last_q .eq. last_s) .and. &
         (switch_pq .eqv. switch_rs)) then
!
         call eri%construct_g_symm_from_L(L_J_pq_p, g_pqrs, alpha, beta, dim_p, dim_q)
!
      else
         call eri%L_pointer_setup_t1(L_J_rs_p, switch_rs, &
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
   end subroutine construct_g_pqrs_t1_t1_eri_tool_c
!
!
   pure subroutine get_eri_t1_packed_mem_t1_eri_tool_c(eri, string, req_pq, &
                                                       dim_p, dim_q, qp)
!!
!!    Get ERI T1 packed mem
!!    Written by Rolf H. Myhre Jun 2020
!!
!!    Adds the memory usage for get ERI that depends on 
!!    dimensions p and q to req_pq
!!    Note that req_pq must be initialized outside
!!    This routine is meant for estimates for the batching manager.
!!    If you are batching over index q, dim_q will typically be 1
!!    and the other dimensions should be full.
!!
!!    Temporary arrays needed for reordering, triggered by optional qp,
!!    are not allocated simultaneously in construct_g, so this routine may
!!    overestimate total memory usage
!!
!!    See get_eri_t1 for additional documentation
!!
      implicit none
!
      class(t1_eri_tool_c), intent(in) :: eri
!
      character(len=2), intent(in) :: string
!
      integer, intent(inout) :: req_pq
!
      integer, intent(in) :: dim_p, dim_q
!
      logical, optional, intent(in) :: qp
!
!     No extra memory needed if integrals in mem
      if(.not. (eri%t1_eri_mem)) then
!
         call eri%get_packed_eri_mem(string, req_pq, dim_p, dim_q, qp)
!
      endif
!
   end subroutine get_eri_t1_packed_mem_t1_eri_tool_c
!
!
   subroutine get_eri_t1_packed_t1_eri_tool_c(eri, string, g_pqpq,              &
                                              first_p, last_p, first_q, last_q, &
                                              alpha, beta, qp)
!!
!!    Get ERI T1 packed
!!    Written by Rolf H. Myhre, Jan 2021
!!
!!    based on routines by Sarai D. Folkestad and Eirik F. Kjønstad
!!
!!    Returns packed eri integrals in the rectangular full packed (RFP) normal upper format
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
!!    Make sure to test both when using this routine.
!!
!!    See: ACM Transactions on Mathematical Software, Vol. 37, No. 2, Article 18
!!         http://doi.acm.org/10.1145/1731022.1731028
!!
!!    string: two character string indicating integral block. Example: "vo"
!!            o: occupied, v: virtual, f: full
!!
!!    q_pqpq: array to contain the integral
!!
!!    first_p, last_p, etc.: first and last index of integrals
!!                           in range determined by string
!!
!!    alpha: scales data added to g_pqrs (default = 1.0)
!!    beta : scales data already in g_pqrs (default = 0.0)
!!
!!    qp:   switches order of p and q, i.e., reorders to g_qpqp
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      character(len=2), intent(in) :: string
!
      complex(dp), intent(inout), dimension(1) :: g_pqpq
!
      integer, optional, intent(in) :: first_p, last_p
      integer, optional, intent(in) :: first_q, last_q
!
      complex(dp), optional, intent(in) :: alpha, beta
!
      logical, optional, intent(in) :: qp
!
      complex(dp) :: alpha_, beta_
!
      integer :: full_first_p, full_last_p, full_first_q, full_last_q, dim_p, dim_q
!
      alpha_ = one_complex
      if(present(alpha)) alpha_ = alpha
!
      beta_ = zero_complex
      if(present(beta)) beta_ = beta
!
      call eri%index_setup(string(1:1), full_first_p, full_last_p, first_p, last_p)
      call eri%index_setup(string(2:2), full_first_q, full_last_q, first_q, last_q)
!
      dim_p = full_last_p - full_first_p + 1
      dim_q = full_last_q - full_first_q + 1
!
      call eri%get_g_pqrs_t1_packed(g_pqpq, full_first_p, full_last_p,               &
                                            full_first_q, full_last_q, dim_p, dim_q, &
                                            alpha_, beta_, qp)
!
   end subroutine get_eri_t1_packed_t1_eri_tool_c
!
!
   subroutine get_g_pqrs_t1_packed_t1_eri_tool_c(eri, g_pqpq, first_p, last_p, first_q, last_q, &
                                                              dim_p, dim_q, alpha, beta, qp)
!!
!!    Get g_pqrs T1 packed
!!    Written by Rolf H. Myhre, Jan 2021
!!
!!    based on routines by Sarai D. Folkestad and Eirik F. Kjønstad
!!
!!    q_pqpq: array to contain the integral
!!    first_p, last_p, etc.: absolute first and last index of integrals
!!    dim_p, dim_q: length of p and q range
!!
!!    Optional reordering logicals (default = .false.):
!!    qp:   switches order of p and q, i.e., g_qpqp
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
      integer, intent(in) :: dim_p, dim_q
!
      complex(dp), intent(inout), dimension(dim_p*dim_q*(dim_p*dim_q+1)/2) :: g_pqpq
!
      complex(dp), intent(in) :: alpha, beta
!
      logical, optional, intent(in) :: qp
!
      if (eri%t1_eri_mem) then
!
         call eri%copy_g_pqrs_to_packed(g_pqpq, eri%g_pqrs_t1, &
                                        first_p, first_q,      &
                                        dim_p, dim_q, alpha, beta, qp)
!
      else
!
         call eri%construct_g_pqrs_t1_packed(g_pqpq, first_p, last_p, first_q, last_q, &
                                             dim_p, dim_q, alpha, beta, qp)
!
      endif
!
!
   end subroutine get_g_pqrs_t1_packed_t1_eri_tool_c
!
!
   subroutine construct_g_pqrs_t1_packed_t1_eri_tool_c(eri, g_pqpq, &
                                                       first_p, last_p, first_q, last_q, &
                                                       dim_p, dim_q, alpha, beta, qp)
!!
!!    Construct g_pqrs T1 packed
!!    Written by Rolf H. Myhre, Jan 2021
!!
!!    q_pqpq: array to contain the integral
!!    first_p, last_p, etc.: absolute first and last index of integrals
!!    dim_p, dim_q: length of p and q range
!!
!!    Optional reordering logicals (default = .false.):
!!    qp:   switches order of p and q, i.e., reorders to g_qpqp
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
      integer, intent(in) :: dim_p, dim_q
!
      complex(dp), intent(out), dimension(dim_p*dim_q*(dim_p*dim_q+1)/2) :: g_pqpq
!
      complex(dp), intent(in) :: alpha, beta
!
      logical, optional, intent(in) :: qp
!
      complex(dp), dimension(:,:,:), pointer, contiguous :: L_J_pq_p
!
      logical :: switch_pq,  pq_alloced
!
      L_J_pq_p => null()
!
      switch_pq = .false.
      if(present(qp)) switch_pq = qp
!
      call eri%L_pointer_setup_t1(L_J_pq_p, switch_pq, &
                                  first_p, last_p, first_q, last_q, pq_alloced)
!
      call eri%construct_g_packed_from_L(L_J_pq_p, g_pqpq, alpha, beta, dim_p, dim_q)
!
      if (pq_alloced) then
         call mem%dealloc(L_J_pq_p, eri%n_J, dim_p, dim_q)
      endif
!
      L_J_pq_p => null()
!
   end subroutine construct_g_pqrs_t1_packed_t1_eri_tool_c
!
!
   pure subroutine get_eri_c1_mem_t1_eri_tool_c(eri, string, &
                                                req_p, req_q, req_r, req_s, req_pq, req_rs, &
                                                dim_p, dim_q, dim_r, dim_s,  &
                                                qp, sr)
!!
!!    Get ERI C1 mem
!!    Written by Rolf H. Myhre Jun 2020
!!
!!    Adds the memory usage for get ERI C1 depending on pq and rs to req_pq and req_rs
!!    Note that req_pq and req_rs must be initialized outside
!!    This routine is meant for estimates for the batching manager.
!!    If you are batching over index r, dim_r will typically be 1
!!    and the other dimensions should be full.
!!
!!    Temporary arrays needed for reordering, triggered by optional qp and sr
!!    are not allocated simultaneously in construct_eri, so this routine may
!!    overestimate total memory usage
!!
!!    See get_eri_c1 for additional documentation
!!
      implicit none
!
      class(t1_eri_tool_c), intent(in) :: eri
!
      character(len=4), intent(in) :: string
!
      integer, intent(inout) :: req_pq, req_rs
      integer, intent(inout) :: req_p, req_q, req_r, req_s
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      logical, optional, intent(in) :: qp, sr
!
      integer :: req_pq_c1, req_rs_c1
      integer :: req_pq_t1, req_rs_t1
!
      req_pq_c1 = eri%n_J*dim_p*dim_q
      req_rs_c1 = eri%n_J*dim_r*dim_s
      req_pq_t1 = 0
      req_rs_t1 = 0
!
      call eri%get_eri_t1_mem(string, req_pq_t1, req_rs_t1, dim_p, dim_q, dim_r, dim_s, qp, sr)
!
      call eri%construct_cholesky_c1_mem(string(1:2), req_p, req_q, req_pq_c1, dim_p, dim_q, qp)
      call eri%construct_cholesky_c1_mem(string(3:4), req_r, req_s, req_rs_c1, dim_r, dim_s, sr)
!
      req_pq = req_pq + max(req_pq_t1, req_pq_c1)
      req_rs = req_rs + max(req_rs_t1, req_rs_c1)
!
   end subroutine get_eri_c1_mem_t1_eri_tool_c
!
!
   subroutine get_eri_c1_t1_eri_tool_c(eri, string, g_pqrs, c_ai,        &
                                       first_p, last_p, first_q, last_q, &
                                       first_r, last_r, first_s, last_s, &
                                       alpha, beta, qp, sr, rspq)
!!
!!    Get ERI C1
!!    Written by Rolf H. Myhre Jun 2020
!!
!!    Returns the integral like tensor [G_pqrs, C_ai]
!!
!!    g_pqrs_c1 = sum_J L_J_pq_t1 L_J_rs_c1 + sum_J L_J_pq_c1 L_J_rs_t1
!!
!!    string: four character string indicating integral block. Example: "vovo"
!!            o: occupied, v: virtual, f: full
!!    q_pqrs: array to contain the integral
!!    c_ai  : T1 like array used in the construction
!!    first_p, last_p, etc.: first and last index of integrals in range determined by string
!!
!!    alpha: scales data added to g_pqrs (default = 1.0)
!!    beta : scales data already in g_pqrs (default = 0.0)
!!
!!    qp:   optional logical, switches order of p and q, i.e., g_qprs, default: false
!!    sr:   optional logical, switches order of r and s, i.e., g_pqsr, default: false
!!    rspq: optional logical, switches pq and rs, i.e., g_rspq, default: false
!!    Note that all reordering switches can be combined freely
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      character(len=4), intent(in) :: string
!
      complex(dp), dimension(eri%n_v, eri%n_o), intent(in) :: c_ai
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
      call eri%construct_g_pqrs_c1(g_pqrs, c_ai, &
                                   full_first_p, full_last_p, &
                                   full_first_q, full_last_q, &
                                   full_first_r, full_last_r, &
                                   full_first_s, full_last_s, &
                                   alpha_, beta_, qp, sr, rspq)
!
   end subroutine get_eri_c1_t1_eri_tool_c
!
!
   subroutine construct_g_pqrs_c1_t1_eri_tool_c(eri, g_pqrs, c_ai,                &
                                                first_p, last_p, first_q, last_q, &
                                                first_r, last_r, first_s, last_s, &
                                                alpha, beta, qp, sr, rspq)
!!
!!    Construct g_pqrs C1
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    based on routines by Alexander C. Paul
!!
!!    Constructs the integral like tensor [g_pqrs, C_ai]
!!
!!    g_pqrs_c1 = sum_J L_J_pq_t1 L_J_rs_c1 + sum_J L_J_pq_c1 L_J_rs_t1
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
      class(t1_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
      integer, intent(in) :: first_r, last_r
      integer, intent(in) :: first_s, last_s
!
      complex(dp), intent(inout), &
                dimension(first_p:last_p,first_q:last_q,first_r:last_r,first_s:last_s) :: g_pqrs
!
      complex(dp), dimension(eri%n_v, eri%n_o), intent(in) :: c_ai
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

      logical :: alloced
!
      complex(dp) :: beta_one
!
      L_J_pq_p => null()
      L_J_rs_p => null()
!
      switch_pq = .false.
      switch_rs = .false.
      if(present(qp)) switch_pq = qp
      if(present(sr)) switch_rs = sr
!
      beta_one = beta
!
      if ((last_p .le. eri%n_o) .and. (last_q .gt. eri%n_o) .and. &
          (last_r .le. eri%n_o) .and. (last_s .gt. eri%n_o)) then
         call output%error_msg('tried to construct ovov_c1')
      endif
!
      dim_p = last_p - first_p + 1
      dim_q = last_q - first_q + 1
      dim_r = last_r - first_r + 1
      dim_s = last_s - first_s + 1
!
      if (.not. ((last_p .le. eri%n_o) .and. (last_q .gt. eri%n_o))) then
!
         call eri%L_pointer_setup_t1(L_J_rs_p, switch_rs, first_r,last_r, first_s,last_s, alloced)
!
         call mem%alloc(L_J_pq_p, eri%n_J, dim_p, dim_q)
         call eri%construct_cholesky_c1(L_J_pq_p, c_ai, first_p,last_p, first_q,last_q, switch_pq)
!
         call eri%construct_g_from_L(L_J_pq_p, L_J_rs_p, g_pqrs, alpha, beta, &
                                     dim_p, dim_q, dim_r, dim_s, rspq)
!
         beta_one = one
!
         call mem%dealloc(L_J_pq_p, eri%n_J, dim_p, dim_q)
         if (alloced) then
            call mem%dealloc(L_J_rs_p, eri%n_J, dim_r, dim_s)
         else
            L_J_rs_p => null()
         endif
      endif
!
      if (.not. ((last_r .le. eri%n_o) .and. (last_s .gt. eri%n_o))) then
!
         call eri%L_pointer_setup_t1(L_J_pq_p, switch_pq, first_p,last_p, first_q,last_q, alloced)
!
         call mem%alloc(L_J_rs_p, eri%n_J, dim_r, dim_s)
         call eri%construct_cholesky_c1(L_J_rs_p, c_ai, first_r,last_r, first_s,last_s, switch_rs)
!
         call eri%construct_g_from_L(L_J_pq_p, L_J_rs_p, g_pqrs, alpha, beta_one, &
                                     dim_p, dim_q, dim_r, dim_s, rspq)
!
         call mem%dealloc(L_J_rs_p, eri%n_J, dim_r, dim_s)
         if (alloced) then
            call mem%dealloc(L_J_pq_p, eri%n_J, dim_p, dim_q)
         else
            L_J_pq_p => null()
         endif
!
      endif
!
   end subroutine construct_g_pqrs_c1_t1_eri_tool_c
!
!
   pure subroutine construct_cholesky_c1_mem_t1_eri_tool_c(eri, string, req_p, req_q, req_pq, &
                                                           dim_p, dim_q, qp)
!!
!!    Construct Cholesky C1 mem
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    Figure out how much memory is needed to construct a "C1 Cholesky vector"
!!
!!    See construct_cholesky_c1 for additional documentation
!!
      implicit none
!
      class(t1_eri_tool_c), intent(in) :: eri
      character(len=2), intent(in)     :: string
      integer, intent(inout)           :: req_p, req_q, req_pq
      integer, intent(in)              :: dim_p, dim_q
      logical, intent(in), optional    :: qp
!
      logical :: switch_pq
!
      switch_pq = .false.
      if(present(qp)) switch_pq = qp
!
      if (string(1:1) .eq. 'o') then
!
         if(.not. eri%cholesky_mem) req_p = req_p + eri%n_J*dim_p*eri%n_v
         if(switch_pq)              req_pq = req_pq + eri%n_J*dim_p*dim_q
!
      elseif (string(1:1) .eq. 'v') then ! zero for ov
!
         if (string(2:2) .eq. 'o') then ! vo
!
            req_p = req_p + eri%n_J*eri%n_o*dim_q
            req_pq = req_pq + eri%n_J*dim_p*dim_q

            if(.not. switch_pq .and. dim_p .lt. eri%n_o) then !Will overestimate if batching over p
               req_q = req_q + eri%n_J*eri%n_o*dim_q
            endif
!
            if(.not. eri%cholesky_mem) then
               req_p = req_p + eri%n_J*eri%n_v*dim_p
            endif
!
            if(switch_pq .and. .not. eri%cholesky_mem) then
               req_q = req_q + eri%n_J*eri%n_o*dim_q
            endif
!
         else ! vv
!
            req_q = req_q + eri%n_J*eri%n_o*dim_q
            if (.not. switch_pq) then
               req_pq = req_pq + eri%n_J*dim_q*dim_p
               if(dim_p .lt. eri%n_o) then !Will overestimate if batching over p
                  req_q = req_q + eri%n_J*dim_q*eri%n_o
               endif
            elseif(.not. eri%cholesky_mem) then
               req_q = req_q + eri%n_J*dim_q*eri%n_o
            endif
!
         endif
      endif
!
   end subroutine construct_cholesky_c1_mem_t1_eri_tool_c
!
!
   subroutine construct_cholesky_c1_t1_eri_tool_c(eri, L_J_pq_c1, c_ai, &
                                                  first_p, last_p, first_q, last_q, qp)
!!
!!    Construct Cholesky C1
!!
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    Constructs the "C1 Cholesky vector"
!!
!!    qp: optional (default = .false) reordering logical, returns L_J_qp_c1 if .true.
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      complex(dp), dimension(eri%n_J, first_p:last_p, first_q:last_q), intent(out) :: L_J_pq_c1
!
      complex(dp), dimension(eri%n_v, eri%n_o), intent(in) :: c_ai
!
      logical, intent(in) :: qp
!
      if (.not. eri%is_pq_in_block(first_p, last_p, first_q, last_q)) then
            call output%error_msg('Tried to construct c1 cholesky not in block')
      endif
!
      if (last_q .le. eri%n_o) then
         if (last_p .le. eri%n_o) then
            call eri%construct_cholesky_oo_c1(L_J_pq_c1, c_ai, &
                                              first_p, last_p, first_q, last_q, qp)
         else
            call eri%construct_cholesky_vo_c1(L_J_pq_c1, c_ai, &
                                              first_p-eri%n_o, last_p-eri%n_o, &
                                              first_q, last_q, qp)
         endif
      else
         if (last_p .gt. eri%n_o) then
            call eri%construct_cholesky_vv_c1(L_J_pq_c1, c_ai, &
                                              first_p-eri%n_o, last_p-eri%n_o, &
                                              first_q-eri%n_o, last_q-eri%n_o, qp)
         else
            call output%error_msg('Tried to construct ov c1 Cholesky')
         endif
      endif
!
   end subroutine construct_cholesky_c1_t1_eri_tool_c
!
!
   subroutine construct_cholesky_oo_c1_t1_eri_tool_c(eri, L_J_ij_c1, c_bj, &
                                                     first_i, last_i, first_j, last_j, qp)
!!
!!    Construct Cholesky oo C1
!!    Written by Rolf. H. Myhre, Jun 2020
!!
!!    based on routines by Alexander C. Paul
!!
!!    Computes
!!
!!       L_J_ij_C1 = sum_b L_J_ib_T1 c_bj
!!
!!    and returns the result in L_J_ij_C1
!!
!!    qp: optional (default = .false) reordering logical, returns L_J_ji_C1 if .true.
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_i, last_i
      integer, intent(in) :: first_j, last_j
!
      complex(dp), dimension(eri%n_J,first_i:last_i,first_j:last_j), &
                   target, intent(out) :: L_J_ij_c1
!
      complex(dp), dimension(eri%n_v, eri%n_o), intent(in) :: c_bj
!
      logical, intent(in) :: qp
!
      complex(dp), dimension(:,:,:), allocatable :: L_J_ib
      complex(dp), dimension(:,:,:), pointer :: L_J_ij_p
!
      integer :: i_length, j_length
!
      L_J_ij_p => null()
      i_length = last_i - first_i + 1
      j_length = last_j - first_j + 1
!
      if(.not. qp) then
         L_J_ij_p => L_J_ij_c1
      else
         call mem%alloc(L_J_ij_p, eri%n_J, i_length, j_length)
      endif
!
      if (eri%cholesky_mem) then
!
!        L_J_ij_c1 = sum_b L_J_ib_t1*c1_bj
         call zgemm('N','N',                                      &
                     i_length*eri%n_J, j_length, eri%n_v,         &
                     one_complex,                                 &
                     eri%L_J_ov_t1(:,first_i,1), eri%n_J*eri%n_o, & ! L_Ji_b
                     c_bj(:, first_j), eri%n_v,                   & ! c_b_j
                     zero_complex,                                &
                     L_J_ij_p,                                    & ! L_Ji_j
                     i_length*eri%n_J)
!
      else
!
!        Read the t1-transformed Cholesky vectors
!
         call mem%alloc(L_J_ib, eri%n_J, i_length, eri%n_v)
         call eri%get_cholesky_t1(L_J_ib, first_i, last_i, eri%n_o + 1, eri%n_mo)
!
!        L_J_ij_c1 = sum_b L_J_ib_t1*c1_bj
         call zgemm('N','N',                              &
                     i_length*eri%n_J, j_length, eri%n_v, &
                     one_complex,                         &
                     L_J_ib, i_length*eri%n_J,            & ! L_Ji_b
                     c_bj(:, first_j), eri%n_v,           & ! c_b_j
                     zero_complex,                        &
                     L_J_ij_p,                            & ! L_Ji_j
                     i_length*eri%n_J)
!
         call mem%dealloc(L_J_ib, eri%n_J, i_length, eri%n_v)
!
      endif
!
      if(qp) then
         call sort_123_to_132(L_J_ij_p, L_J_ij_c1, eri%n_J, i_length, j_length)
         call mem%dealloc(L_J_ij_p, eri%n_J, i_length, j_length)
      endif
!
   end subroutine construct_cholesky_oo_c1_t1_eri_tool_c
!
!
   subroutine construct_cholesky_vo_c1_t1_eri_tool_c(eri, L_J_ai_c1, c_ai, &
                                                     first_a, last_a, first_i, last_i, qp)
!!
!!    Construct Cholesky vo C1
!!    Written by Rolf. H. Myhre, Jun 2020
!!
!!    based on routines by Alexander C. Paul
!!
!!    Computes
!!
!!       L_J_ai_C1 = sum_b L_J_ab_T1 c_bi - sum_j c_aj L_J_ji_T1
!!
!!    and returns the result in L_J_ai_c1
!!
!!    qp: optional (default = .false) reordering logical, returns L_J_ia_C1 if .true.
!!
      use array_utilities, only: zero_array_complex
!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      integer, intent(in) :: first_i, last_i
      integer, intent(in) :: first_a, last_a
!
      complex(dp), dimension(eri%n_J,first_a:last_a,first_i:last_i), &
                   target, intent(out) :: L_J_ai_c1
!
      complex(dp), dimension(eri%n_v, eri%n_o), intent(in) :: c_ai
!
      logical, intent(in) :: qp
!
      complex(dp), dimension(:,:,:), allocatable         :: L_J_ij, L_J_ab
      complex(dp), dimension(:,:,:), allocatable, target :: L_J_ia
      complex(dp), dimension(:,:,:), pointer             :: L_J_ia_p
!
      integer :: length_a, length_i
!
      length_i = last_i - first_i + 1
      length_a = last_a - first_a + 1
!
      call mem%alloc(L_J_ij, eri%n_J, length_i, eri%n_o)
!
      if(.not. qp) then
         call mem%alloc(L_J_ia, eri%n_J, length_i, max(length_a, eri%n_o))
         L_J_ia_p => L_J_ia
      else
         if(.not. eri%cholesky_mem) call mem%alloc(L_J_ia, eri%n_J, eri%n_o, length_i)
         L_J_ia_p => L_J_ai_c1
      endif
!
      if (eri%cholesky_mem) then
         call sort_123_to_132(eri%L_J_oo_t1(:,:,first_i:last_i), L_J_ij, eri%n_J,eri%n_o,length_i)
      else
         call eri%get_cholesky_t1(L_J_ia, 1, eri%n_o, first_i, last_i)
         call sort_123_to_132(L_J_ia, L_J_ij, eri%n_J, eri%n_o, length_i)
      endif
!
      call zgemm('N','T',                              &
                  eri%n_J*length_i, length_a, eri%n_o, &
                  -one_complex,                        &
                  L_J_ij, eri%n_J*length_i,            & !L_J_ij
                  c_ai(first_a:, 1), eri%n_v,          & !C_a_j
                  zero_complex,                        &
                  L_J_ia_p, eri%n_J*length_i)            !L_J_ia
!
      call mem%dealloc(L_J_ij, eri%n_J, length_i, eri%n_o)
!
      if(.not. qp) then
         call sort_123_to_132(L_J_ia_p, L_J_ai_c1, eri%n_J, length_i, length_a)
         call mem%dealloc(L_J_ia, eri%n_J, length_i, max(length_a, eri%n_o))
         L_J_ia_p => L_J_ai_c1
      else
         if(.not. eri%cholesky_mem) call mem%dealloc(L_J_ia, eri%n_J, eri%n_o, length_i)
         call mem%alloc(L_J_ia, eri%n_J, length_i, length_a)
         call zero_array_complex(L_J_ia, eri%n_J*length_i*length_a)
         L_J_ia_p => L_J_ia
      endif
!
      if (eri%cholesky_mem) then
!
         call zgemm('N', 'N',                                    &
                    eri%n_J*length_a, length_i, eri%n_v,         &
                    one_complex,                                 &
                    eri%L_J_vv_t1(:,first_a,1), eri%n_J*eri%n_v, & ! L_Ja_b
                    c_ai(:, first_i), eri%n_v,                   & ! c_b_i
                    one_complex,                                 &
                    L_J_ia_p, eri%n_J*length_a)                    ! L_J_ai
!
      else
!
         call mem%alloc(L_J_ab, eri%n_J, length_a, eri%n_v)
         call eri%get_cholesky_t1(L_J_ab, eri%n_o+first_a, eri%n_o+last_a, eri%n_o + 1, eri%n_mo)
!
         call zgemm('N', 'N',                            &
                    eri%n_J*length_a, length_i, eri%n_v, &
                    one_complex,                         &
                    L_J_ab, eri%n_J*length_a,            & ! L_Ja_b
                    c_ai(:, first_i), eri%n_v,           & ! c_b_i
                    one_complex,                         &
                    L_J_ia_p, eri%n_J*length_a)            ! L_Ja_i
!
         call mem%dealloc(L_J_ab, eri%n_J, length_a, eri%n_v)
!
      endif
!
      if(qp) then
         call add_132_to_123(one_complex, L_J_ia_p, L_J_ai_c1, eri%n_J, length_i, length_a)
         call mem%dealloc(L_J_ia, eri%n_J, length_i, length_a)
      endif
!
   end subroutine construct_cholesky_vo_c1_t1_eri_tool_c
!
!
   subroutine construct_cholesky_vv_c1_t1_eri_tool_c(eri, L_J_ab_c1, c_ai, &
                                                     first_a, last_a, first_b, last_b, qp)
!!
!!    Construct Cholesky vv C1
!!    Written by Rolf. H. Myhre, Jun 2020
!!
!!    based on routines by Alexander C. Paul
!!
!!    Computes
!!
!!       L_ab_J_c1= - sum_i c_ai L_ib_J_T1 ,
!!
!!    and returns the result in L_J_ab_c1
!!
!!    qp: optional (default = .false) reordering logical, returns L_J_ba_c1 if .true.
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
!
      complex(dp), dimension(eri%n_v, eri%n_o), intent(in) :: c_ai
!
      integer, intent(in) :: first_a, last_a
      integer, intent(in) :: first_b, last_b
!
      complex(dp), dimension(eri%n_J,first_a:last_a,first_b:last_b), target, intent(out) :: L_J_ab_c1
!
      logical, intent(in) :: qp
!
      complex(dp), dimension(:,:,:), allocatable         :: L_J_bi
      complex(dp), dimension(:,:,:), allocatable, target :: L_J_ba
      complex(dp), dimension(:,:,:), pointer             :: L_J_ba_p
!
      integer :: length_b, length_a
!
      length_a = last_a - first_a + 1
      length_b = last_b - first_b + 1
!
      call mem%alloc(L_J_bi, eri%n_J, length_b, eri%n_o)
!
      if(.not. qp) then
         call mem%alloc(L_J_ba, eri%n_J, length_b, max(length_a,eri%n_o))
         L_J_ba_p => L_J_ba
      else
         if(.not. eri%cholesky_mem) call mem%alloc(L_J_ba, eri%n_J, length_b, eri%n_o)
         L_J_ba_p => L_J_ab_c1
      endif
!
      if (eri%cholesky_mem) then
!
         call sort_123_to_132(eri%L_J_ov_t1(:,:,first_b:last_b), L_J_bi, eri%n_J,eri%n_o,length_b)
!
      else
!
         call eri%get_cholesky_t1(L_J_ba, 1, eri%n_o, eri%n_o + first_b, eri%n_o + last_b)
         call sort_123_to_132(L_J_ba, L_J_bi, eri%n_J, eri%n_o, length_b)
!
      endif
!
!     Calculate and add c1-transformed L_ab_J_c1 = - sum_i c_ai L_ib_J
!
      call zgemm('N','T',                              &
                  length_b*eri%n_J, length_a, eri%n_o, &
                  -one_complex,                        &
                  L_J_bi, eri%n_J*length_b,            & !L_J_bi
                  c_ai(first_a:, 1), eri%n_v,          & !C_a_i
                  zero_complex,                        &
                  L_J_ba_p, eri%n_J*length_b)            !L_J_ba
!
      if(.not. qp) then
         call sort_123_to_132(L_J_ba_p, L_J_ab_c1, eri%n_J, length_b, length_a)
      endif
!
      L_J_ba_p => null()
      call mem%dealloc(L_J_bi, eri%n_J, length_b, eri%n_o)
      if(.not. qp) then
         call mem%dealloc(L_J_ba, eri%n_J, length_b, max(length_a,eri%n_o))
      else
         if(.not. eri%cholesky_mem) call mem%dealloc(L_J_ba, eri%n_J, length_b, eri%n_o)
      endif
!
   end subroutine construct_cholesky_vv_c1_t1_eri_tool_c
!
!
   subroutine set_pointer_t1_t1_eri_tool_c(eri, point, last_p, last_q, first_q)
!!
!!    set pointer T1
!!    written by Rolf H. Myhre, Jun 2020
!!
!!    Set pointer to Cholesky vector block
!!    Assumes full index p
!!
      implicit none
!
      class(t1_eri_tool_c), intent(in), target :: eri
      complex(dp), dimension(:,:,:), pointer, contiguous, intent(out) :: point
      integer, intent(in) :: last_p, last_q, first_q
!
      if(associated(point)) call output%error_msg('Pointer associated in set pointer t1')
      if(.not. eri%cholesky_mem) call output%error_msg('Vectors not in mem for set poiner t1')
!
      if (last_q .le. eri%n_o) then
         if (last_p .le. eri%n_o) then
            point => eri%L_J_oo_t1(:,:,first_q:last_q)
         else
            point => eri%L_J_vo_t1(:,:,first_q:last_q)
         endif
      else
         if (last_p .le. eri%n_o) then
            point => eri%L_J_ov_t1(:,:,first_q-eri%n_o:last_q-eri%n_o)
         else
            point => eri%L_J_vv_t1(:,:,first_q-eri%n_o:last_q-eri%n_o)
         endif
      endif
!
   end subroutine set_pointer_t1_t1_eri_tool_c
!
!
   subroutine L_pointer_setup_t1_t1_eri_tool_c(eri, point, switch, &
                                               first_p, last_p, first_q, last_q, alloced)
!!
!!    pointer setup T1
!!    written by Rolf H. Myhre, Jun 2020
!!
!!    Do not use this routine!
!!
!!    Sets the pointer point to either point to a stored Cholesky vector
!!    or allocates the pointer and fills it with a vector in appropriate order
!!
!!    switch: logical, point will point to L_J_pq instead of L_J_qp if .true.
!!
!!    Deallocing point if alloced = .false. will take you on exciting adventures
!!
      implicit none
!
      class(t1_eri_tool_c), intent(inout) :: eri
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
            call eri%set_pointer_t1(point, last_p, last_q, first_q)
         else
            call eri%get_cholesky_t1(point, first_p, last_p, first_q, last_q)
         endif
!
      else !Reordering
!
!        Set temp to vectors if possible, else read from file
!
         if(vec_in_mem) then
            call eri%set_pointer_t1(temp, last_p, last_q, first_q)
         else
            call eri%get_cholesky_t1(temp, first_p, last_p, first_q, last_q)
         endif
!
!        Sort into point
!
         call sort_123_to_132(temp, point, eri%n_J, last_p-first_p+1, last_q-first_q+1)
!
!        Unset temp if in mem, else dealloc
!
         if (vec_in_mem) then
            temp => null()
         else
            call mem%dealloc(temp, eri%n_J, last_p-first_p+1, last_q-first_q+1)
         endif
!
      endif
!
   end subroutine L_pointer_setup_t1_t1_eri_tool_c
!
!
   subroutine cleanup_t1_eri_tool_c(eri)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2019
!!
      implicit none
!
      class(t1_eri_tool_c) :: eri
!
      if (eri%cholesky_mem) then
!
!        Deallocate Cholesky arrays
!
         call mem%dealloc(eri%L_J_oo_t1, eri%n_J, eri%n_o, eri%n_o)
         call mem%dealloc(eri%L_J_ov_t1, eri%n_J, eri%n_o, eri%n_v)
         call mem%dealloc(eri%L_J_vo_t1, eri%n_J, eri%n_v, eri%n_o)
         call mem%dealloc(eri%L_J_vv_t1, eri%n_J, eri%n_v, eri%n_v)
!
         call mem%dealloc(eri%L_J_oo_mo, eri%n_J, eri%n_o, eri%n_o)
         call mem%dealloc(eri%L_J_ov_mo, eri%n_J, eri%n_o, eri%n_v)
         call mem%dealloc(eri%L_J_vo_mo, eri%n_J, eri%n_v, eri%n_o)
         call mem%dealloc(eri%L_J_vv_mo, eri%n_J, eri%n_v, eri%n_v)
!
      endif
!
!     Deallocate electron repulsion integrals
!
      if (eri%mo_eri_mem) then
         call mem%dealloc(eri%g_pqrs_mo, eri%n_mo, eri%n_mo, &
                                         eri%n_mo, eri%n_mo)
      endif
      if (eri%t1_eri_mem) then
         call mem%dealloc(eri%g_pqrs_t1, eri%n_mo, eri%n_mo, &
                                         eri%n_mo, eri%n_mo)
      endif
!
      eri%cholesky_mem = .false.
      eri%mo_eri_mem = .false.
      eri%t1_eri_mem = .false.
!
   end subroutine cleanup_t1_eri_tool_c
!
!
end module t1_eri_tool_c_class
