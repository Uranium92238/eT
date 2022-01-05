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
submodule (uhf_class) set_get_uhf
!
!!
!!    Set get submodule
!!
!!    Gathers routines that set and get the UHF type-bound variables.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine set_ao_fock_uhf(wf, F)
!!
!!    Set AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the AO Fock.
!!
!
      use reordering, only: squareup
!
      implicit none
!
      class(uhf) :: wf
!
      real(dp), dimension(wf%ao%n*(wf%ao%n + 1)/2, wf%n_densities), intent(in) :: F ! Packed
!
      real(dp), dimension(:), allocatable :: F_sigma
!
      call mem%alloc(F_sigma, wf%ao%n*(wf%ao%n + 1)/2)
!
!     Alpha Fock
!
      call dcopy(wf%ao%n*(wf%ao%n + 1)/2, F, 1, F_sigma, 1)
      call squareup(F_sigma, wf%ao_fock_a, wf%ao%n)
!
!     Beta Fock
!
      call dcopy(wf%ao%n*(wf%ao%n + 1)/2, F(1, 2), 1, F_sigma, 1)
      call squareup(F_sigma, wf%ao_fock_b, wf%ao%n)
!
      call mem%dealloc(F_sigma, wf%ao%n*(wf%ao%n + 1)/2)
!
   end subroutine set_ao_fock_uhf
!
!
   module subroutine get_ao_fock_uhf(wf, F)
!!
!!    Get AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Returns the AO Fock
!!
!
      use reordering, only: packin
!
      implicit none
!
      class(uhf), intent(in) :: wf
!
      real(dp), dimension(wf%ao%n*(wf%ao%n+1)/2, wf%n_densities), intent(inout) :: F ! Packed
!
      real(dp), dimension(:), allocatable :: F_sigma
!
      call mem%alloc(F_sigma, wf%ao%n*(wf%ao%n + 1)/2)
!
!     Alpha Fock
!
      call packin(F_sigma, wf%ao_fock_a, wf%ao%n)
      call dcopy(wf%ao%n*(wf%ao%n + 1)/2, F_sigma, 1, F, 1)
!
!     Beta Fock
!
      call packin(F_sigma, wf%ao_fock_b, wf%ao%n)
      call dcopy(wf%ao%n*(wf%ao%n + 1)/2, F_sigma, 1, F(1, 2), 1)
!
      call mem%dealloc(F_sigma, wf%ao%n*(wf%ao%n + 1)/2)
!
   end subroutine get_ao_fock_uhf
!
!
   module subroutine get_ao_density_sq_uhf(wf, D)
!!
!!    Get AO density squared
!!    Written by Eirik F. Kjønstad, Nov 2018
!!
!!    Returns the unpacked AO density matrix D
!!    (or density matrices in descendants, see overwriting routines)
!!
      implicit none
!
      class(uhf), intent(in) :: wf
!
      real(dp), dimension(wf%ao%n**2, wf%n_densities), intent(inout) :: D
!
      call dcopy(wf%ao%n**2, wf%ao_density_a, 1, D, 1)
      call dcopy(wf%ao%n**2, wf%ao_density_b, 1, D(1, 2), 1)
!
   end subroutine get_ao_density_sq_uhf
!
!
end submodule set_get_uhf
