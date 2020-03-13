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
   module subroutine get_ovov_ccs(wf, g_iajb, first_i, last_i, first_a, last_a, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_iajb
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_b, last_b
!
   end subroutine get_ovov_ccs
!
!
   module subroutine get_oooo_ccs(wf, g_ijkl, first_i, last_i, first_j, last_j, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ijkl
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_k, last_k
      integer, optional, intent(in) :: first_l, last_l
!
   end subroutine get_oooo_ccs
!
!
   module subroutine get_ooov_ccs(wf, g_ijka, first_i, last_i, first_j, last_j, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ijka
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_k, last_k
      integer, optional, intent(in) :: first_a, last_a
!
   end subroutine get_ooov_ccs
!
!
   module subroutine get_oovo_ccs(wf, g_ijak, first_i, last_i, first_j, last_j, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ijak
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_k, last_k
      integer, optional, intent(in) :: first_a, last_a
!
   end subroutine get_oovo_ccs
!
!
   module subroutine get_ovoo_ccs(wf, g_iajk, first_i, last_i, first_a, last_a, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_iajk
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_k, last_k
      integer, optional, intent(in) :: first_a, last_a
!
   end subroutine get_ovoo_ccs
!
!
   module subroutine get_vooo_ccs(wf, g_aijk, first_a, last_a, first_i, last_i, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_aijk
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_k, last_k
      integer, optional, intent(in) :: first_a, last_a
!
   end subroutine get_vooo_ccs
!
!
   module subroutine get_vvoo_ccs(wf, g_abij, first_a, last_a, first_b, last_b, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_abij
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
   end subroutine get_vvoo_ccs
!
!
   module subroutine get_vovo_ccs(wf, g_aibj, first_a, last_a, first_i, last_i, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_aibj
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
   end subroutine get_vovo_ccs
!
!
   module subroutine get_voov_ccs(wf, g_aijb, first_a, last_a, first_i, last_i, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_aijb
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
   end subroutine get_voov_ccs
!
!
   module subroutine get_ovvo_ccs(wf, g_iabj, first_i, last_i, first_a, last_a, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_iabj
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
   end subroutine get_ovvo_ccs
!
!
   module subroutine get_oovv_ccs(wf, g_ijab, first_i, last_i, first_j, last_j, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ijab
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
   end subroutine get_oovv_ccs
!
!
   module subroutine get_vvvo_ccs(wf, g_abci, first_a, last_a, first_b, last_b, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_abci
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_c, last_c
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
   end subroutine get_vvvo_ccs
!
!
   module subroutine get_vvov_ccs(wf, g_abic, first_a, last_a, first_b, last_b, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_abic
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_c, last_c
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
   end subroutine get_vvov_ccs
!
!
   module subroutine get_vovv_ccs(wf, g_aibc, first_a, last_a, first_i, last_i, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_aibc
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_c, last_c
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
   end subroutine get_vovv_ccs
!
!
   module subroutine get_ovvv_ccs(wf, g_iabc, first_i, last_i, first_a, last_a, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_iabc
      integer, optional, intent(in) :: first_i, last_i
      integer, optional, intent(in) :: first_c, last_c
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
   end subroutine get_ovvv_ccs
!
!
   module subroutine get_vvvv_ccs(wf, g_abcd, first_a, last_a, first_b, last_b, &
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
      class(ccs), intent(inout) :: wf
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_abcd
      integer, optional, intent(in) :: first_d, last_d
      integer, optional, intent(in) :: first_c, last_c
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
   end subroutine get_vvvv_ccs
!
!
   module subroutine get_g_pqrs_required_ccs(wf, req_l, req_r, dim_p, dim_q, dim_r, dim_s)
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
      integer, intent(in)  :: dim_p, dim_q, dim_r, dim_s
      integer, intent(out) :: req_l, req_r
!
   end subroutine get_g_pqrs_required_ccs
