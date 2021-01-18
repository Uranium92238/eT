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
submodule (hf_class) set_get_hf
!
!!
!!    Set get submodule
!!
!!    Gathers routines that set and get the HF type-bound variables.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine get_ao_density_sq_hf(wf, D)
!!
!!    Get AO density squared
!!    Written by Eirik F. Kjønstad, Nov 2018
!!
!!    Returns the unpacked AO density matrix D
!!    (or density matrices in descendants, see overwriting routines)
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao**2,wf%n_densities), intent(inout) :: D
!
      call dcopy(wf%n_ao**2, wf%ao_density, 1, D, 1)
!
   end subroutine get_ao_density_sq_hf
!
!
   module subroutine set_ao_density_hf(wf, D)
!!
!!    Set AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the AO density from input
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:) :: D ! Packed
!
      call squareup(D, wf%ao_density, wf%n_ao)
!
   end subroutine set_ao_density_hf
!
!
   module subroutine set_ao_fock_hf(wf, F)
!!
!!    Set AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the AO Fock from input
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao*(wf%n_ao + 1)/2, wf%n_densities), intent(in) :: F ! Packed
!
      call squareup(F(:,1), wf%ao_fock, wf%n_ao)
!
   end subroutine set_ao_fock_hf
!
!
   module subroutine get_ao_fock_hf(wf, F)
!!
!!    Set AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the AO Fock from input
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao*(wf%n_ao+1)/2, wf%n_densities), intent(inout) :: F ! Packed
!
      call packin(F(:,1), wf%ao_fock, wf%n_ao)
!
   end subroutine get_ao_fock_hf
!
!
   module subroutine get_ao_density_hf(wf, D)
!!
!!    Get AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Packs the AO density into D.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(:,:) :: D
!
      call packin(D(:,1), wf%ao_density, wf%n_ao)
!
   end subroutine get_ao_density_hf
!
!
end submodule set_get_hf