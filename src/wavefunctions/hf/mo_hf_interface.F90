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
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp) :: max_gradient
!
   end function get_max_roothan_hall_mo_gradient_hf
!
!
!
   module subroutine get_roothan_hall_mo_gradient_hf(wf, G)
!!
!!    Get MO Roothan-Hall gradient
!!    Written by Sarai D. Folkestad, 2018
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o), intent(out) :: G
!
   end subroutine get_roothan_hall_mo_gradient_hf
!
!
   module subroutine get_mo_fock_hf(wf, F)
!!
!!    Get Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(:,:), intent(inout) :: F
!
   end subroutine get_mo_fock_hf
!
!
   module subroutine set_mo_fock_hf(wf, F)
!!
!!    Set Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the MO Fock from input
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
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
      implicit none
!
      class(hf) :: wf
!
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
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
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
!!    Prepare for Roothan-Hall MO
!!    Written by Sarai D. Folkestad, Oct 2019
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine prepare_for_roothan_hall_mo_hf
!