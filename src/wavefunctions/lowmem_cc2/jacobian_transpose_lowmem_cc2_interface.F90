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
   module subroutine effective_jacobian_transpose_transformation_lowmem_cc2(wf, omega, b, cvs)
!!
!!    Effective Jacobian transpose transformation
!!    Written by Sarai Dery Folkestad, Jun 2019
!!
      implicit none
!
      class(lowmem_cc2) :: wf
!
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: b
!
      logical, intent(in) :: cvs
!
   end subroutine effective_jacobian_transpose_transformation_lowmem_cc2
!
!
   module subroutine jacobian_transpose_ccs_b1_lowmem_cc2(wf, sigma_ai, b_ai)
!!
!!    Jacobian transpose B1 (CCS)
!!    Written by Sarai D. Folkestad, Jun 2019
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
   end subroutine jacobian_transpose_ccs_b1_lowmem_cc2
!
!
   module subroutine jacobian_transpose_cc2_a1_lowmem_cc2(wf, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Jacobian transpose A1 (CC2)
!!    Written by Sarai D. Folkestad, Jun 2019
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine jacobian_transpose_cc2_a1_lowmem_cc2
!
!
   module subroutine jacobian_transpose_cc2_b1_lowmem_cc2(wf, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Jacobian transpose B1 (CC2)
!!    Written by Sarai D. Folkestad, Jun 2019
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine jacobian_transpose_cc2_b1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_a1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 A1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_transpose_cc2_a1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_b1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 B1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_transpose_cc2_b1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_c1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 C1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_transpose_cc2_c1_lowmem_cc2
!
!
!
   module subroutine effective_jacobian_transpose_cc2_d1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 D1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
!
   end subroutine effective_jacobian_transpose_cc2_d1_lowmem_cc2
!
   module subroutine effective_jacobian_transpose_cc2_e1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 E1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v   
!
   end subroutine effective_jacobian_transpose_cc2_e1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_f1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 F1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_transpose_cc2_f1_lowmem_cc2
!