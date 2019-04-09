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
   module subroutine effective_jacobian_transpose_transformation_cc3(wf, omega, c)
!!
!!    Jacobian transpose transformation (CC3)
!!    Alexander Paul and Rolf H. Myhre, March 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
   end subroutine effective_jacobian_transpose_transformation_cc3
!
!
   module subroutine jacobian_transpose_cc3_A_cc3(wf, omega, c_ai, c_abij, sigma_ai, sigma_abij)
!!
!!    Terms of the transpose of the  CC3 Jacobi matrix
!!    Alex C. Paul and Rolf H. Myhre, March 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: sigma_abij
!
   end subroutine jacobian_transpose_cc3_A_cc3
!
!