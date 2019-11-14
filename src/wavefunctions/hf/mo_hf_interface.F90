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
   module function get_max_roothan_hall_mo_gradient_hf(wf) result(max_gradient)
!!
!!    Get max of Roothan-Hall gradient
!!    Written by Sarai D. Folkestad, 2018
!!
!!    Constructs the Roothan-Hall gradient,
!!
!!       E^(1) = - 4 F_ai,
!!
!!    and returns the maximum absolute value of E^(1).
!!
      implicit none
!
      class(hf), intent(in) :: wf
      real(dp) :: max_gradient
   end function get_max_roothan_hall_mo_gradient_hf
!
!
   module subroutine get_roothan_hall_mo_gradient_hf(wf, G)
!!
!!    Get packed Roothan-Hall gradient
!!    Written by Sarai D. Folkestad, 2018
!!
!!    Constructs the Roothan-Hall gradient,
!!
!!       G_ai = E^(1)_ai = - 4 F_ai
!!
      implicit none
!
      class(hf), intent(in) :: wf
      real(dp), dimension(wf%n_v,wf%n_o), intent(out) :: G
!
   end subroutine get_roothan_hall_mo_gradient_hf
!
!
   module subroutine get_mo_fock_hf(wf, F)
!!
!!    Get MO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Gets the MO Fock
!!
      implicit none
!
      class(hf), intent(in) :: wf
      real(dp), dimension(:,:), intent(inout) :: F
!
   end subroutine get_mo_fock_hf
!
!
   module subroutine set_mo_fock_hf(wf, F)
!!
!!    Set MO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the MO Fock from input
!!
      implicit none
!
      class(hf), intent(in) :: wf
      real(dp), dimension(:, :), intent(inout) :: F
!
   end subroutine set_mo_fock_hf
!
!
   module subroutine roothan_hall_update_orbitals_mo_hf(wf)
!!
!!    Roothan-Hall update of orbitals
!!    Written by Linda Goletto and Sarai D. Folkestad, 2019
!!
!!    Update orbitals after a Roothan-Hall step in the MO basis
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine roothan_hall_update_orbitals_mo_hf
!
!
   module subroutine do_roothan_hall_mo_hf(wf, F, e)
!!
!!    Do Roothan-Hall
!!    Written by Linda Goletto and Sarai D. Folkestad, 2019
!!
!!    Solves the equation 
!!
!!       F W = W e, 
!!
!!    where W (called mo_basis_update) is the 
!!    transformation matrix from old to new MO basis
!!
!!       C_new = C_old W
!!
!!    'F' : Fock matrix in MO basis
!!
!!    'e' : orbital energies 
!!
      implicit none
!
      class(hf) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in)    :: F
      real(dp), dimension(wf%n_mo),          intent(inout) :: e
!
   end subroutine do_roothan_hall_mo_hf
!
!
   module subroutine update_fock_and_energy_mo_hf(wf, h_wx, prev_ao_density)
!!
!!    Update Fock and energy
!!    Written by Linda Goletto and Sarai D. Folkestad, 2019
!!
!!    Constructs Fock matrix in the AO basis
!!    and calculates the energy.
!!
!!    Transforms the active Fock matrix to the new MO basis.
!!
      implicit none
!
      class(hf) :: wf
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
      real(dp), dimension(wf%n_ao**2, wf%n_densities), intent(in), optional :: prev_ao_density
!
   end subroutine update_fock_and_energy_mo_hf
!
!
   module subroutine initialize_W_mo_update_hf(wf)
!!
!!    Initialize W MO update
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad 
!!    and Linda Goletto, Jan 2019
!!
!!    Modified by Ida-Marie Hoyvik, Oct 2019
!!
!!    Initializes the eigenvectors W 
!!    for Roothan-Hall in the mo basis (FW = We)
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine initialize_W_mo_update_hf
!
!
   module subroutine destruct_W_mo_update_hf(wf)
!!
!!    Destruct W MO update 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad 
!!    and Linda Goletto, Jan 2019
!!
!!    Destructs the eigenvectors W 
!!    for Roothan-Hall in the mo basis (FW = We)
!!
!!    Modified by Ida-Marie Hoyvik, Oct 2019
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine destruct_W_mo_update_hf
!
!
   module subroutine prepare_for_roothan_hall_mo_hf(wf)
!!
!!    Construct idempotent density and prepare for Roothan-Hall
!!    Written by Sarai D. Folkestad, Oct 2019
!!
!!    Constructs the ao Fock matrix and
!!    performs a Roothan-Hall step to get the
!!    initial idempotent density. 
!!
!!    The routine also prints the number of 
!!    electrons and the energy of the initial guess.
!!
!!    These steps are performed in 
!!
!!       prepare_for_roothan_hall_hf,
!!
!!    this routine also allocates the MO 
!!    quantities mo_fock and W_mo_update and
!!    initializew W_mo_update to the identity matrix
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine prepare_for_roothan_hall_mo_hf
