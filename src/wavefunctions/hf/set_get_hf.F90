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
      real(dp), dimension(wf%ao%n**2,wf%n_densities), intent(inout) :: D
!
      call dcopy(wf%ao%n**2, wf%ao_density, 1, D, 1)
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
      call squareup(D, wf%ao_density, wf%ao%n)
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
      real(dp), dimension(wf%ao%n*(wf%ao%n + 1)/2, wf%n_densities), intent(in) :: F ! Packed
!
      call squareup(F(:,1), wf%ao_fock, wf%ao%n)
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
      real(dp), dimension(wf%ao%n*(wf%ao%n+1)/2, wf%n_densities), intent(inout) :: F ! Packed
!
      call packin(F(:,1), wf%ao_fock, wf%ao%n)
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
      call packin(D(:,1), wf%ao_density, wf%ao%n)
!
   end subroutine get_ao_density_hf
!
!
   module subroutine get_ao_h_hf(wf, h)
!!
!!    Get AO h
!!    Written by Tor S. Haugland, Oct 2021
!!
!!    Get the Hamiltonian one-electron integrals h in
!!    the AO basis.
!!
      implicit none
!
      class(hf), intent(in) :: wf
      real(dp), dimension(wf%ao%n, wf%ao%n) :: h
!
      call dcopy(wf%ao%n**2, wf%ao%h, 1, h, 1)
!
   end subroutine get_ao_h_hf
!
!
   module subroutine get_ao_g_hf(wf, g, A, B, C, D, precision_, skip)
!!
!!    Get AO g
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020
!!    Adapted to HF by Tor S. Haugland, Oct 2021
!!
!!    Get the electron repulsion integrals g for the
!!    shell quartet (A, B, C, D).
!!
!!    Two optional arguments:
!!
!!       precision_:  (intent in) Double precision real corresponding to the Libint precision
!!                   'epsilon' to use when calculating the integral. Does not guarantee a precision
!!                   to the given value and should therefore be selected conservatively.
!!
!!       skip:       (intent out) If present, this integer will be 1 if Libint decided not to
!!                   calculate the integral; it will be zero otherwise. If it is present, g will
!!                   not be zeroed out if Libint decides not to calculate g. Thus, only pass 'skip'
!!                   to the routine if you wish to avoid zeroing out elements that are negligible.
!!
!!
      implicit none
!
      class(hf), intent(in) :: wf
      integer :: A, B, C, D
      real(dp), dimension(wf%ao%shells(A)%length * wf%ao%shells(B)%length * &
                          wf%ao%shells(C)%length * wf%ao%shells(D)%length), intent(out) :: g
      real(dp), optional, intent(in) :: precision_
      integer, optional, intent(out) :: skip
!
      call wf%ao%get_eri(g, A, B, C, D, precision_, skip)
!
   end subroutine get_ao_g_hf
!
!
end submodule set_get_hf