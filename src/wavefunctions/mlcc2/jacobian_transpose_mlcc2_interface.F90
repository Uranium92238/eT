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
   module subroutine jacobian_transpose_transformation_mlcc2(wf, b)
!!
!!    Jacobian transpose transformation (MLCC2)
!!    Adapted by Sarai D. Folkestad, Jul 2019
!!
!!    Calculates the transpose Jacobian transformation, i.e., the transformation
!!    by the transpose of the Jacobian matrix
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | R >.
!!
!!    The transformation is performed as sigma^T = c^T A, where c is the vector
!!    sent to the routine. On exit, the vector c is equal to sigma (the transformed
!!    vector).
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: b
!
   end subroutine jacobian_transpose_transformation_mlcc2
!
!
   module subroutine jacobian_transpose_cc2_a1_mlcc2(wf, sigma_ai, c_ai, n_cc2_o, n_cc2_v, &
                                                      first_o, first_v)
!!
!!    Jacobian transpose MLCC2 A1
!!    Written by Sarai D. Folkestad, Apr 2019
!!
!!    Adapted from jacobian_transpose_cc2.F90
!!    written by Sarai D. Folkestad and Alexander C. Paul, Feb 2019
!!
!!    Calculates the A1 term,
!!
!!       A1: sum_bcjk u^bc_jk (c_bj L_iakc - c_ak g_jbic - c_ci g_jbka)
!!
!!    and adds it to sigma_ai.
!!
!!    Index restrictions:
!!
!!       b, c, j, k : CC2 orbitals
!!
!!       a, i : unrestricted
!!
      implicit none
!
      class(mlcc2) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)     :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_o, first_v
!
   end subroutine jacobian_transpose_cc2_a1_mlcc2
!
!
   module subroutine jacobian_transpose_cc2_b1_mlcc2(wf, sigma_ai, c_aibj, n_cc2_o, n_cc2_v, &
                                                      first_o, first_v, last_o, last_v)
!!
!!    Jacobian transpose MLCC2 B1
!!    Written by Sarai D. Folkestad, Apr 2019
!!
!!    Adapted from jacobian_transpose_cc2.F90
!!    written by Sarai D. Folkestad and Alexander C. Paul, Feb 2019
!!
!!    Calculates the B1 term,
!!
!!       B1: sum_bjc c_bjci g_bjca - sum_bjk c_akbj g_bjik
!!
!!    and adds it to sigma_ai.
!!
!!    Index restrictions:
!!
!!       Term 1: 
!!
!!          b, j, c, i : CC2 orbitals
!!
!!          a : unretricted
!!
!!       Term 2: 
!!
!!          a, k, b, j : CC2 orbitals
!!
!!          i : unretricted
!!
      implicit none
!
      class(mlcc2) :: wf
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_o, first_v, last_o, last_v
      real(dp), dimension(n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
!
   end subroutine jacobian_transpose_cc2_b1_mlcc2
!
!
   module subroutine jacobian_transpose_cc2_a2_mlcc2(wf, sigma_aibj, c_ai, n_cc2_o, n_cc2_v, &
                                                      first_o, first_v, last_o, last_v)
!!
!!    Jacobian transpose MLCC2 A2
!!    Written by Sarai D. Folkestad, Apr 2019
!!
!!    Adapted from jacobian_transpose_cc2.F90
!!    written by Sarai D. Folkestad and Alexander C. Paul, Feb 2019
!!
!!    Calculates the A2 term,
!!
!!       A2:  2F_jb c_ai - F_ib c_aj - L_ikjb c_ak + L_cajb c_ci
!!
!!    and adds it to sigma_aibj.
!!
!!    Term 4 is calculated in batches of index c.
!!
!!    Index restrictions:
!!
!!       a, i, b, j : CC2 orbitals
!!
!!       Term 3:
!!
!!          k : unrestricted 
!!
!!       Term 4:
!!
!!          c : unrestricted 
!!
      implicit none
!
      class(mlcc2) :: wf
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_o, first_v, last_o, last_v
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                          :: c_ai
      real(dp), dimension(n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o), intent(inout)   :: sigma_aibj
!
   end subroutine jacobian_transpose_cc2_a2_mlcc2
!
!
   module subroutine jacobian_transpose_cc2_b2_mlcc2(wf, sigma_aibj, c_aibj)
!!
!!    Jacobian transpose CC2 B2
!!    Written by Sarai D. Folkestad, Apr 2019
!!
!!    Adapted from jacobian_transpose_cc2.F90
!!    written by Sarai D. Folkestad and Alexander C. Paul, Feb 2019
!!
!!    Calculates the B2 term,
!!
!!       B2: ε_aibj c_aibj
!!
!!    and adds it to sigma_aibj.
!!
!!    Index restrictions:
!!
!!       a, i, b, j : CC2 orbitals 
!!
      implicit none
!
      class(mlcc2) :: wf
      real(dp), dimension(wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o), intent(in)      :: c_aibj
      real(dp), dimension(wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o), intent(inout)   :: sigma_aibj
!
   end subroutine jacobian_transpose_cc2_b2_mlcc2
