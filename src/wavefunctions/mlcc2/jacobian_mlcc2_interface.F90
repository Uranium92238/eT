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
   module subroutine prepare_for_jacobian_mlcc2(wf)
!!
!!    Prepare for Jacobian
!!    Adapted by Sarai D. Folkestad, 2019
!!
      implicit none 
!
      class(mlcc2), intent(inout) :: wf 
!
   end subroutine prepare_for_jacobian_mlcc2
!
!
   module subroutine jacobian_transformation_mlcc2(wf, c)
!!
!!    Jacobian transformation (mlcc2)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, Jan 2019
!!
!!    Directs the transformation by the cc2 Jacobi matrix,
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | nu >,
!!
!!    where the basis employed for the brackets is biorthonormal.
!!    The transformation is rho = A c, i.e.,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck
!!                  + 1/2 sum_ckdl A_mu,ckdl c_ckdl (1 + delta_ck,dl).
!!
!!    On exit, c is overwritten by rho. That is, c_a_i = rho_a_i,
!!    and c_aibj = rho_aibj.
!!
      implicit none
!
      class(mlcc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
   end subroutine jacobian_transformation_mlcc2
!
!
   module subroutine jacobian_cc2_a1_mlcc2(wf, rho_ai, c_ai, n_cc2_o, n_cc2_v, &
                                          first_cc2_o, first_cc2_v)
!!
!!    Jacobian CC2 A1
!!    Adapted by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlcc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)   :: rho_ai
!
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v
!
   end subroutine jacobian_cc2_a1_mlcc2
!
!
   module subroutine jacobian_cc2_b1_mlcc2(wf, rho_ai, c_aibj, n_cc2_o, n_cc2_v, &
                                          first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v)
!!
!!    Jacobian CC2 B1
!!    Adapted by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v
!
      real(dp), dimension(n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o), intent(in)  :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)                     :: rho_ai   
!
   end subroutine jacobian_cc2_b1_mlcc2
!
!
   module subroutine jacobian_cc2_a2_mlcc2(wf, rho_aibj, c_ai, n_cc2_o, n_cc2_v, &
                                          first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v)
!!
!!    Jacobian CC2 A2
!!    Adapted by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                                  :: c_ai
      real(dp), dimension(wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o), intent(out) :: rho_aibj   
!
   end subroutine jacobian_cc2_a2_mlcc2
!
!
   module subroutine jacobian_cc2_b2_mlcc2(wf, rho_aibj, c_aibj)
!!
!!    Jacobian MLCC2 B2
!!    Adapted by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o), intent(inout)   :: c_aibj
      real(dp), dimension(wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o), intent(inout)  :: rho_aibj   
!
   end subroutine jacobian_cc2_b2_mlcc2
!
!
   module subroutine save_jacobian_a1_intermediates_mlcc2(wf, n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v)
!!
!!    Save jacobian a1 intermediates
!!    Adapted by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v
!
   end subroutine save_jacobian_a1_intermediates_mlcc2
