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
submodule (ccs_class) tei_ccs_complex
!
!!
!!    Two-electron integrals submodule (CCS)
!!    Set up by Andreas Skeidsvoll, Aug 2019
!!
!!    Submodule containing routines that can be used to construct t1-transformed two-electron integrals.
!!
!
      implicit none
!
!
contains
!
!
   module subroutine get_ovov_ccs_complex(wf, g_iajb, first_i, last_i, first_a, last_a, &
                                         first_j, last_j, first_b, last_b)
!!
!!    Get ovov
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_iajb
!
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i
      integer :: local_first_a, local_last_a
      integer :: local_first_j, local_last_j
      integer :: local_first_b, local_last_b
!
      logical :: index_restrictions
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_a) .and. present(last_a) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_b) .and. present(last_b)) then
!
         index_restrictions = .true.
!
         local_first_i = first_i
         local_first_a = first_a
         local_first_j = first_j
         local_first_b = first_b
!
         local_last_i = last_i
         local_last_a = last_a
         local_last_j = last_j
         local_last_b = last_b
!
      else
!
         index_restrictions = .false.
!
         local_first_i = 1
         local_first_a = 1
         local_first_j = 1
         local_first_b = 1
!
         local_last_i = wf%n_o
         local_last_a = wf%n_v
         local_last_j = wf%n_o
         local_last_b = wf%n_v
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_iajb, &
                          local_first_i, local_last_i, &
                          wf%n_o + local_first_a, wf%n_o + local_last_a, &
                          local_first_j, local_last_j, &
                          wf%n_o + local_first_b, wf%n_o + local_last_b)
!
   end subroutine get_ovov_ccs_complex
!
!
   module subroutine get_oooo_ccs_complex(wf, g_ijkl, first_i, last_i, first_j, last_j, &
                                         first_k, last_k, first_l, last_l)
!!
!!    Get oooo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ijkl
!
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_k, last_k
      integer, optional, intent(in) :: first_l, last_l
!
      integer :: local_first_i, local_last_i
      integer :: local_first_j, local_last_j
      integer :: local_first_k, local_last_k
      integer :: local_first_l, local_last_l
!
      logical :: index_restrictions
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_k) .and. present(last_k) .and. &
          present(first_l) .and. present(last_l)) then
!
         index_restrictions = .true.
!
         local_first_i = first_i
         local_first_j = first_j
         local_first_k = first_k
         local_first_l = first_l
!
         local_last_i = last_i
         local_last_j = last_j
         local_last_k = last_k
         local_last_l = last_l
!
      else
!
         index_restrictions = .false.
!
         local_first_i = 1
         local_first_j = 1
         local_first_k = 1
         local_first_l = 1
!
         local_last_i = wf%n_o
         local_last_j = wf%n_o
         local_last_k = wf%n_o
         local_last_l = wf%n_o
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_ijkl, &
                          local_first_i, local_last_i, &
                          local_first_j, local_last_j, &
                          local_first_k, local_last_k, &
                          local_first_l, local_last_l)
!
   end subroutine get_oooo_ccs_complex
!
!
   module subroutine get_ooov_ccs_complex(wf, g_ijka, first_i, last_i, first_j, last_j, &
                                         first_k, last_k, first_a, last_a)
!!
!!    Get ooov
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ijka
!
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_k, last_k
      integer, optional, intent(in) :: first_a, last_a
!
      integer :: local_first_i, local_last_i
      integer :: local_first_j, local_last_j
      integer :: local_first_k, local_last_k
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_k) .and. present(last_k) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i
         local_first_j = first_j
         local_first_k = first_k
         local_first_a = first_a
!
         local_last_i = last_i
         local_last_j = last_j
         local_last_k = last_k
         local_last_a = last_a
!
      else
!
         local_first_i = 1
         local_first_j = 1
         local_first_k = 1
         local_first_a = 1
!
         local_last_i = wf%n_o
         local_last_j = wf%n_o
         local_last_k = wf%n_o
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_ijka, &
                          local_first_i, local_last_i, &
                          local_first_j, local_last_j, &
                          local_first_k, local_last_k, &
                          wf%n_o + local_first_a, wf%n_o + local_last_a)
!
   end subroutine get_ooov_ccs_complex
!
!
   module subroutine get_oovo_ccs_complex(wf, g_ijak, first_i, last_i, first_j, last_j, &
                                         first_a, last_a, first_k, last_k)
!!
!!    Get oovo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ijak
!
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_k, last_k
      integer, optional, intent(in) :: first_a, last_a
!
      integer :: local_first_i, local_last_i
      integer :: local_first_j, local_last_j
      integer :: local_first_k, local_last_k
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_k) .and. present(last_k) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i
         local_first_j = first_j
         local_first_k = first_k
         local_first_a = first_a
!
         local_last_i = last_i
         local_last_j = last_j
         local_last_k = last_k
         local_last_a = last_a
!
      else
!
         local_first_i = 1
         local_first_j = 1
         local_first_k = 1
         local_first_a = 1
!
         local_last_i = wf%n_o
         local_last_j = wf%n_o
         local_last_k = wf%n_o
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_ijak, &
                           local_first_i, local_last_i, &
                           local_first_j, local_last_j, &
                           wf%n_o + local_first_a, wf%n_o + local_last_a, &
                           local_first_k, local_last_k)
!
   end subroutine get_oovo_ccs_complex
!
!
   module subroutine get_ovoo_ccs_complex(wf, g_iajk, first_i, last_i, first_a, last_a, &
                                         first_j, last_j, first_k, last_k)
!!
!!    Get ovoo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_iajk
!
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_k, last_k
      integer, optional, intent(in) :: first_a, last_a
!
      integer :: local_first_i, local_last_i
      integer :: local_first_j, local_last_j
      integer :: local_first_k, local_last_k
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_k) .and. present(last_k) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i
         local_first_j = first_j
         local_first_k = first_k
         local_first_a = first_a
!
         local_last_i = last_i
         local_last_j = last_j
         local_last_k = last_k
         local_last_a = last_a
!
      else
!
         local_first_i = 1
         local_first_j = 1
         local_first_k = 1
         local_first_a = 1
!
         local_last_i = wf%n_o
         local_last_j = wf%n_o
         local_last_k = wf%n_o
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_iajk, &
                           local_first_i, local_last_i, &
                           wf%n_o + local_first_a, wf%n_o + local_last_a, &
                           local_first_j, local_last_j, &
                           local_first_k, local_last_k)
!
   end subroutine get_ovoo_ccs_complex
!
!
   module subroutine get_vooo_ccs_complex(wf, g_aijk, first_a, last_a, first_i, last_i, &
                                         first_j, last_j, first_k, last_k)
!!
!!    Get vooo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_aijk
!
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_k, last_k
      integer, optional, intent(in) :: first_a, last_a
!
      integer :: local_first_i, local_last_i
      integer :: local_first_j, local_last_j
      integer :: local_first_k, local_last_k
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_k) .and. present(last_k) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i
         local_first_j = first_j
         local_first_k = first_k
         local_first_a = first_a
!
         local_last_i = last_i
         local_last_j = last_j
         local_last_k = last_k
         local_last_a = last_a
!
      else
!
         local_first_i = 1
         local_first_j = 1
         local_first_k = 1
         local_first_a = 1
!
         local_last_i = wf%n_o
         local_last_j = wf%n_o
         local_last_k = wf%n_o
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_aijk, &
                           wf%n_o + local_first_a, wf%n_o + local_last_a, &
                           local_first_i, local_last_i, &
                           local_first_j, local_last_j, &
                           local_first_k, local_last_k)
!
   end subroutine get_vooo_ccs_complex
!
!
   module subroutine get_vvoo_ccs_complex(wf, g_abij, first_a, last_a, first_b, last_b, &
                                         first_i, last_i, first_j, last_j)
!!
!!    Get vvoo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_abij
!
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i
      integer :: local_first_j, local_last_j
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i
         local_first_j = first_j
         local_first_b = first_b
         local_first_a = first_a
!
         local_last_i = last_i
         local_last_j = last_j
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1
         local_first_j = 1
         local_first_b = 1
         local_first_a = 1
!
         local_last_i = wf%n_o
         local_last_j = wf%n_o
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_abij, &
                           wf%n_o + local_first_a, wf%n_o + local_last_a, &
                           wf%n_o + local_first_b, wf%n_o + local_last_b, &
                           local_first_i, local_last_i, &
                           local_first_j, local_last_j)
!
   end subroutine get_vvoo_ccs_complex
!
!
   module subroutine get_vovo_ccs_complex(wf, g_aibj, first_a, last_a, first_i, last_i, &
                                         first_b, last_b, first_j, last_j)
!!
!!    Get vovo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_aibj
!
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i
      integer :: local_first_j, local_last_j
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      logical :: index_restrictions
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         index_restrictions = .true.
!
         local_first_i = first_i
         local_first_j = first_j
         local_first_b = first_b
         local_first_a = first_a
!
         local_last_i = last_i
         local_last_j = last_j
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         index_restrictions = .false.
!
         local_first_i = 1
         local_first_j = 1
         local_first_b = 1
         local_first_a = 1
!
         local_last_i = wf%n_o
         local_last_j = wf%n_o
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_aibj, &
                           wf%n_o + local_first_a, wf%n_o + local_last_a, &
                           local_first_i, local_last_i, &
                           wf%n_o + local_first_b, wf%n_o + local_last_b, &
                           local_first_j, local_last_j)
!
   end subroutine get_vovo_ccs_complex
!
!
   module subroutine get_voov_ccs_complex(wf, g_aijb, first_a, last_a, first_i, last_i, &
                                         first_j, last_j, first_b, last_b)
!!
!!    Get voov
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_aijb
!
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i
      integer :: local_first_j, local_last_j
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i
         local_first_j = first_j
         local_first_b = first_b
         local_first_a = first_a
!
         local_last_i = last_i
         local_last_j = last_j
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1
         local_first_j = 1
         local_first_b = 1
         local_first_a = 1
!
         local_last_i = wf%n_o
         local_last_j = wf%n_o
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_aijb, &
                           wf%n_o + local_first_a, wf%n_o + local_last_a, &
                           local_first_i, local_last_i, &
                           local_first_j, local_last_j, &
                           wf%n_o + local_first_b, wf%n_o + local_last_b)
!
   end subroutine get_voov_ccs_complex
!
!
   module subroutine get_ovvo_ccs_complex(wf, g_iabj, first_i, last_i, first_a, last_a, &
                                         first_b, last_b, first_j, last_j)
!!
!!    Get ovvo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_iabj
!
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i
      integer :: local_first_j, local_last_j
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i
         local_first_j = first_j
         local_first_b = first_b
         local_first_a = first_a
!
         local_last_i = last_i
         local_last_j = last_j
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1
         local_first_j = 1
         local_first_b = 1
         local_first_a = 1
!
         local_last_i = wf%n_o
         local_last_j = wf%n_o
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_iabj, &
                           local_first_i, local_last_i, &
                           wf%n_o + local_first_a, wf%n_o + local_last_a, &
                           wf%n_o + local_first_b, wf%n_o + local_last_b, &
                           local_first_j, local_last_j)
!
   end subroutine get_ovvo_ccs_complex
!
!
   module subroutine get_oovv_ccs_complex(wf, g_ijab, first_i, last_i, first_j, last_j, &
                                         first_a, last_a, first_b, last_b)
!!
!!    Get oovv
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ijab
!
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i
      integer :: local_first_j, local_last_j
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i
         local_first_j = first_j
         local_first_b = first_b
         local_first_a = first_a
!
         local_last_i = last_i
         local_last_j = last_j
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1
         local_first_j = 1
         local_first_b = 1
         local_first_a = 1
!
         local_last_i = wf%n_o
         local_last_j = wf%n_o
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_ijab, &
                           local_first_i, &
                           local_last_i, &
                           local_first_j, &
                           local_last_j, &
                           wf%n_o + local_first_a, &
                           wf%n_o + local_last_a, &
                           wf%n_o + local_first_b, &
                           wf%n_o + local_last_b)
!
   end subroutine get_oovv_ccs_complex
!
!
   module subroutine get_vvvo_ccs_complex(wf, g_abci, first_a, last_a, first_b, last_b, &
                                       first_c, last_c, first_i, last_i)
!!
!!    Get vvvo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_abci
!
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_c, last_c
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i
      integer :: local_first_c, local_last_c
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_c) .and. present(last_c) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i
         local_first_c = first_c
         local_first_b = first_b
         local_first_a = first_a
!
         local_last_i = last_i
         local_last_c = last_c
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1
         local_first_c = 1
         local_first_b = 1
         local_first_a = 1
!
         local_last_i = wf%n_o
         local_last_c = wf%n_v
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_abci, &
                           wf%n_o + local_first_a, &
                           wf%n_o + local_last_a, &
                           wf%n_o + local_first_b, &
                           wf%n_o + local_last_b, &
                           wf%n_o + local_first_c, &
                           wf%n_o + local_last_c,&
                           local_first_i, &
                           local_last_i)
!
   end subroutine get_vvvo_ccs_complex
!
!
   module subroutine get_vvov_ccs_complex(wf, g_abic, first_a, last_a, first_b, last_b, &
                                         first_i, last_i, first_c, last_c)
!!
!!    Get vvov
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_abic
!
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_c, last_c
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i
      integer :: local_first_c, local_last_c
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_c) .and. present(last_c) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i
         local_first_c = first_c
         local_first_b = first_b
         local_first_a = first_a
!
         local_last_i = last_i
         local_last_c = last_c
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1
         local_first_c = 1
         local_first_b = 1
         local_first_a = 1
!
         local_last_i = wf%n_o
         local_last_c = wf%n_v
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_abic, &
                           wf%n_o + local_first_a, &
                           wf%n_o + local_last_a, &
                           wf%n_o + local_first_b, &
                           wf%n_o + local_last_b, &
                           local_first_i, &
                           local_last_i, &
                           wf%n_o + local_first_c, &
                           wf%n_o + local_last_c)
!
   end subroutine get_vvov_ccs_complex
!
!
   module subroutine get_vovv_ccs_complex(wf, g_aibc, first_a, last_a, first_i, last_i, &
                                         first_b, last_b, first_c, last_c)
!!
!!    Get vovv
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_aibc
!
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_c, last_c
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i
      integer :: local_first_c, local_last_c
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_c) .and. present(last_c) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i
         local_first_c = first_c
         local_first_b = first_b
         local_first_a = first_a
!
         local_last_i = last_i
         local_last_c = last_c
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1
         local_first_c = 1
         local_first_b = 1
         local_first_a = 1
!
         local_last_i = wf%n_o
         local_last_c = wf%n_v
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_aibc, &
                           wf%n_o + local_first_a, &
                           wf%n_o + local_last_a, &
                           local_first_i, &
                           local_last_i, &
                           wf%n_o + local_first_b, &
                           wf%n_o + local_last_b, &
                           wf%n_o + local_first_c, &
                           wf%n_o + local_last_c)
!
   end subroutine get_vovv_ccs_complex
!
!
   module subroutine get_ovvv_ccs_complex(wf, g_iabc, first_i, last_i, first_a, last_a, &
                                         first_b, last_b, first_c, last_c)
!!
!!    Get ovvv
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_iabc
!
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_c, last_c
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i
      integer :: local_first_c, local_last_c
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_c) .and. present(last_c) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i
         local_first_c = first_c
         local_first_b = first_b
         local_first_a = first_a
!
         local_last_i = last_i
         local_last_c = last_c
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1
         local_first_c = 1
         local_first_b = 1
         local_first_a = 1
!
         local_last_i = wf%n_o
         local_last_c = wf%n_v
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_iabc, &
                           local_first_i, &
                           local_last_i, &
                           wf%n_o + local_first_a, &
                           wf%n_o + local_last_a, &
                           wf%n_o + local_first_b, &
                           wf%n_o + local_last_b, &
                           wf%n_o + local_first_c, &
                           wf%n_o + local_last_c)
!
   end subroutine get_ovvv_ccs_complex
!
!
   module subroutine get_vvvv_ccs_complex(wf, g_abcd, first_a, last_a, first_b, last_b, &
                                         first_c, last_c, first_d, last_d)
!!
!!    Get vvvv
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions
!!    are provided, the routines assume that the full integral should be returned.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_abcd
!
      integer, optional, intent(in) :: first_d, last_d
      integer, optional, intent(in) :: first_c, last_c
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_d, local_last_d
      integer :: local_first_c, local_last_c
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      logical :: index_restrictions
!
      if (present(first_d) .and. present(last_d) .and. &
          present(first_c) .and. present(last_c) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         index_restrictions = .true.
!
         local_first_d = first_d
         local_first_c = first_c
         local_first_b = first_b
         local_first_a = first_a
!
         local_last_d = last_d
         local_last_c = last_c
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         index_restrictions = .false.
!
         local_first_d = 1
         local_first_c = 1
         local_first_b = 1
         local_first_a = 1
!
         local_last_d = wf%n_v
         local_last_c = wf%n_v
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%get_g_pqrs_t1_complex(g_abcd, &
                           wf%n_o + local_first_a, &
                           wf%n_o + local_last_a, &
                           wf%n_o + local_first_b, &
                           wf%n_o + local_last_b, &
                           wf%n_o + local_first_c, &
                           wf%n_o + local_last_c, &
                           wf%n_o + local_first_d, &
                           wf%n_o + local_last_d)
!
   end subroutine get_vvvv_ccs_complex
!
!
   module subroutine get_g_pqrs_required_ccs_complex(wf, req_l, req_r, dim_p, dim_q, dim_r, dim_s)
!!
!!    Get memory required to construct g_pqrs
!!    Written by Rolf H. Myhre, April 2019
!!
!!    Simple routine calculate an integral block with provided dimensions.
!!    req_l and req_r are the memory required by construct_g_pqrs to allocate 
!!    the Cholesky vectors
!!
!!    req_l = n_J*dim_p*dim_q, left Cholesky vector
!!    req_r = n_J*dim_r*dim_s, right Cholesky vector
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      integer, intent(in)  :: dim_p, dim_q, dim_r, dim_s
      integer, intent(out) :: req_l, req_r
!
      req_l = wf%integrals%n_J*dim_p*dim_q
      req_r = wf%integrals%n_J*dim_r*dim_s
!
   end subroutine get_g_pqrs_required_ccs_complex
!
!
end submodule tei_ccs_complex
