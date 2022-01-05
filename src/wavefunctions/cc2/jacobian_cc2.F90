!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
submodule (cc2_class) jacobian_cc2
!
!!
!!    Jacobian submodule
!!
!!    Routines for the linear transform of trial
!!    vectors by the Jacobian matrix
!!
!!    ρ_i = A * c_i,
!!
!!    where
!!
!!    A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_jacobian_cc2(wf)
!!
!!    Prepare for jacobian
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      type(timings), allocatable :: timer
!
      timer = timings('Prepare for Jacobian CC2', pl='normal')
      call timer%turn_on()
!
      call wf%initialize_t2()
      call wf%construct_t2()
!
      call wf%doubles%prepare_for_jacobian()
!
      call timer%turn_off()
!
   end subroutine prepare_for_jacobian_cc2
!
!
   module subroutine jacobian_doubles_b2_cc2(wf, rho_aibj, c_aibj)
!!
!!    Jacobian doubles B2 CC2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2019 and 2021
!!
!!    rho_aibj =+ c_aibj * (eps_a - eps_i)
!!
!!    Note that a symmetrization is needed after this term.
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)    :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: rho_aibj
!
      integer :: i, j, a, b
!
      type(timings) :: timer
!
      timer = timings('Jacobian doubles B2 CC2', pl='verbose')
      call timer%turn_on()
!
!$omp parallel do private(a, i, b, j)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  rho_aibj(a,i,b,j) = rho_aibj(a,i,b,j) + c_aibj(a,i,b,j)*&
                                          (- wf%orbital_energies(i) &
                                           + wf%orbital_energies(wf%n_o + a))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call timer%turn_off()
!
   end subroutine jacobian_doubles_b2_cc2
!
!
end submodule jacobian_cc2
