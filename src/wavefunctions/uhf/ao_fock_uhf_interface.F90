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
   module subroutine construct_mo_fock_uhf(wf)
!!
!!    Construct MO Fock
!!    Written by Eirik F. Kjønstad, Mar 2019
!!
!!    Give notice to user that it does not yet exist.
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
   end subroutine construct_mo_fock_uhf
!
!
   module subroutine update_fock_and_energy_cumulative_uhf(wf, prev_ao_density)
!!
!!    Update Fock and energy cumulatively
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    This routine guides the construction of the Fock matrix (or matrices for
!!    unrestricted wavefunctions) from the current AO density (or densities).
!!    It is called by the solver and is overwritten for unrestricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
      real(dp), dimension(wf%ao%n**2, wf%n_densities), intent(in) :: prev_ao_density
!
   end subroutine update_fock_and_energy_cumulative_uhf
!
!
   module subroutine update_fock_and_energy_non_cumulative_uhf(wf)
!!
!!    Update Fock and energy
!!    Written by Eirik F. Kjønstad, Sep 2018
!!    Modified for QM/MM by Tommaso Giovannini, July 2019 
!!
!!    This routine guides the construction of the Fock matrix (or matrices for
!!    unrestricted wavefunctions) from the current AO density (or densities).
!!    It is called by the solver and is overwritten for unrestricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine update_fock_and_energy_non_cumulative_uhf
!
!
   module subroutine update_fock_and_energy_uhf(wf, prev_ao_density)
!!
!!    Update Fock and energy wrapper
!!    Written by Tommaso Giovannini, July 2019
!!
!!    Wrapper for cumulative or non-cumulative subroutines
!!    depending on the path
!!
      implicit none
!
      class(uhf) :: wf
      real(dp), dimension(wf%ao%n**2, wf%n_densities), intent(in), optional :: prev_ao_density
!
   end subroutine update_fock_and_energy_uhf
!
!
   module subroutine construct_ao_spin_fock_uhf(wf, D, D_sigma, sigma, &
                     h_wx, cumulative)
!!
!!    Construct AO spin Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    The routine computes the alpha or beta Fock matrix, depending
!!    on the value of the spin 'sigma' (='alpha' or 'beta'):
!!
!!       F_αβ^alpha = h_αβ + sum_γδ g_αβγδ D_γδ - sum_γδ g_αδγβ D_γδ^alpha
!!       F_αβ^beta  = h_αβ + sum_γδ g_αβγδ D_γδ - sum_γδ g_αδγβ D_γδ^beta
!!
!!    Here the superscript refers to the spin function, while the subscripts
!!    are AO indices. In contrast to the restricted routine, this one does
!!    not calculate the energy - a separate call is required to get the
!!    unrestricted Hartree-Fock energy.
!!
      implicit none
!
      class(uhf) :: wf
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D_sigma
      character(len=*), intent(in) :: sigma
      logical, intent(in), optional :: cumulative
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: h_wx
!
   end subroutine construct_ao_spin_fock_uhf
