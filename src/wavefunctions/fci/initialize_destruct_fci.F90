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
submodule (fci_class) initialize_destruct
!
!!
!! Initialize destruct submodule
!!
!! Gathers routines that initialize and destruct the FCI type-bound variables.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine initialize_strings(wf)
!!
!!    Initialize strings
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      call mem%alloc(wf%alpha_strings, wf%n_alpha_strings)
      call mem%alloc(wf%beta_strings, wf%n_beta_strings)
!
   end subroutine initialize_strings
!
!
   module subroutine destruct_strings(wf)
!!
!!    Destruct strings
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      call mem%dealloc(wf%alpha_strings, wf%n_alpha_strings)
      call mem%dealloc(wf%beta_strings, wf%n_beta_strings)
!
   end subroutine destruct_strings
!
!
   module subroutine initialize_excitation_maps(wf)
!!
!!    Initialize excitation maps
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      call mem%alloc(wf%excitation_maps_alpha, wf%n_alpha_strings, wf%n_alpha_excitations, 4)
      call mem%alloc(wf%excitation_maps_beta, wf%n_beta_strings, wf%n_beta_excitations, 4)
!
   end subroutine initialize_excitation_maps
!
!
   module subroutine destruct_excitation_maps(wf)
!!
!!    Destruct excitation maps
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      call mem%dealloc(wf%excitation_maps_alpha, wf%n_alpha_strings, wf%n_alpha_excitations, 4)
      call mem%dealloc(wf%excitation_maps_beta, wf%n_beta_strings, wf%n_beta_excitations, 4)
!
   end subroutine destruct_excitation_maps
!
!
   module subroutine initialize_hamiltonian_integrals(wf)
!!
!!    Initialize Hamiltonian integrals
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      call mem%alloc(wf%h_pq, wf%n_mo, wf%n_mo)
      call mem%alloc(wf%g_pqrs, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
      call mem%alloc(wf%effective_2e_hamiltonian, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_hamiltonian_integrals
!
!
   module subroutine destruct_hamiltonian_integrals(wf)
!!
!!    Destruct Hamiltonian integrals
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      call mem%dealloc(wf%h_pq, wf%n_mo, wf%n_mo)
      call mem%dealloc(wf%g_pqrs, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
      call mem%dealloc(wf%effective_2e_hamiltonian, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_hamiltonian_integrals
!
!
   module subroutine initialize_ci_coefficients(wf)
!!
!!    Initialize fci vectors
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      call mem%alloc(wf%ci_coefficients, wf%n_alpha_strings, wf%n_beta_strings, wf%n_states)
!
   end subroutine initialize_ci_coefficients
!
!
   module subroutine destruct_ci_coefficients(wf)
!!
!!    Destruct fci vectors
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      call mem%dealloc(wf%ci_coefficients, wf%n_alpha_strings, wf%n_beta_strings, wf%n_states)
!
   end subroutine destruct_ci_coefficients
!
!
   module subroutine initialize_energies(wf)
!!
!!    Initialize fci energies
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      call mem%alloc(wf%energies, wf%n_states)
!
   end subroutine initialize_energies
!
!
   module subroutine destruct_energies(wf)
!!
!!    Destruct fci energies
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      call mem%dealloc(wf%energies, wf%n_states)
!
   end subroutine destruct_energies
!
!
   module subroutine initialize_gs_density_fci(wf)
!!
!!    Initialize gs density
!!    Written by Alexander C. Paul, May 2022
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      call mem%alloc(wf%density, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_gs_density_fci
!
!
   module subroutine destruct_gs_density_fci(wf)
!!
!!    Destruct gs density
!!    Written by Alexander C. Paul, May 2022
!!
      implicit none
!
      class(fci), intent(inout) :: wf
!
      if (allocated(wf%density)) call mem%dealloc(wf%density, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_gs_density_fci
!
!
end submodule initialize_destruct
