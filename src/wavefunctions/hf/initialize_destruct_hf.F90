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
submodule (hf_class) initialize_destruct_hf
!
!!
!!    Initialize destruct submodule
!!
!!    Gathers routines that initialize and destruct the HF type-bound variables.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine initialize_ao_density_hf(wf)
!!
!!    Initialize AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%ao_density)) call mem%alloc(wf%ao_density, wf%ao%n, wf%ao%n)
      if (.not. allocated(wf%previous_ao_density)) &
         call mem%alloc(wf%previous_ao_density, wf%ao%n, wf%ao%n, wf%n_densities)
!
   end subroutine initialize_ao_density_hf
!
!
   module subroutine destruct_ao_density_hf(wf)
!!
!!    Destruct AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%ao_density)) call mem%dealloc(wf%ao_density, wf%ao%n, wf%ao%n)
      if (allocated(wf%previous_ao_density)) &
         call mem%dealloc(wf%previous_ao_density, wf%ao%n, wf%ao%n, wf%n_densities)
!
   end subroutine destruct_ao_density_hf
!
!
   module subroutine initialize_orbitals_hf(wf)
!!
!!    Initialize orbitals
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the arrays associated with the orbital
!!    coefficients. In spin-unrestricted hfs, this
!!    will include alpha and beta coefficients, though these
!!    are the same and therefore redundant in restricted
!!    hfs.
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
!
      wf%orbital_coefficients = zero
      wf%orbital_energies     = zero
!
   end subroutine initialize_orbitals_hf
!
!
   module subroutine initialize_density_hf(wf)
!!
!!    Initialize density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the AO density (or densities).
!!    In spin-unrestricted hfs, this alpha and beta densities,
!!    though these are the same and therefore redundant in restricted
!!    hfs.
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%initialize_ao_density()
!
      wf%ao_density = zero
!
   end subroutine initialize_density_hf
!
!
   module subroutine initialize_fock_hf(wf)
!!
!!    Initialize Fock
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the AO Fock matrix (or matrices).
!!    In spin-unrestricted hfs, this alpha and beta matrices,
!!    though these are the same and therefore redundant in restricted
!!    hfs.
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%initialize_ao_fock()
!
      wf%ao_fock = zero
!
   end subroutine initialize_fock_hf
!
!
   module subroutine destruct_fock_hf(wf)
!!
!!    Destruct Fock
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Destructs the AO Fock matrix (or matrices).
!!    In spin-unrestricted hfs, this alpha and beta matrices,
!!    though these are the same and therefore redundant in restricted
!!    hfs.
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%destruct_ao_fock()
!
   end subroutine destruct_fock_hf
!
!
   module subroutine initialize_orbital_coefficients_frozen_hf_hf(wf)
!!
!!    Initialize orbital coefficients frozen core
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%orbital_coefficients_frozen_hf)) &
            call mem%alloc(wf%orbital_coefficients_frozen_hf, wf%ao%n, wf%n_frozen_hf_o)
!
   end subroutine initialize_orbital_coefficients_frozen_hf_hf
!
!
   module subroutine destruct_orbital_coefficients_frozen_hf_hf(wf)
!!
!!    Destruct orbital coefficients frozen core
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%orbital_coefficients_frozen_hf)) &
            call mem%dealloc(wf%orbital_coefficients_frozen_hf, wf%ao%n, wf%n_frozen_hf_o)
!
   end subroutine destruct_orbital_coefficients_frozen_hf_hf
!
!
   module subroutine initialize_orbital_coefficients_fc_hf(wf)
!!
!!    Initialize orbital coefficients frozen core
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%orbital_coefficients_fc)) &
            call mem%alloc(wf%orbital_coefficients_fc, wf%ao%n, wf%n_frozen_core_orbitals)
!
   end subroutine initialize_orbital_coefficients_fc_hf
!
!
   module subroutine destruct_orbital_coefficients_fc_hf(wf)
!!
!!    Destruct orbital coefficients frozen core
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%orbital_coefficients_fc)) &
            call mem%dealloc(wf%orbital_coefficients_fc, wf%ao%n, wf%n_frozen_core_orbitals)
!
   end subroutine destruct_orbital_coefficients_fc_hf
!
!
   module subroutine initialize_frozen_CCT_hf(wf)
!!
!!    Initialize frozen CC^T
!!    Written by Sarai D. Folkestad, Jan 2020
!!
!!
      implicit none
!
      class(hf) :: wf
!
      call mem%alloc(wf%frozen_CCT, wf%ao%n, wf%ao%n)
!
   end subroutine initialize_frozen_CCT_hf
!
!
   module subroutine destruct_frozen_CCT_hf(wf)
!!
!!    Destruct frozen CC^T
!!    Written by Sarai D. Folkestad, Jan 2020
!!
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%frozen_CCT)) call mem%dealloc(wf%frozen_CCT, wf%ao%n, wf%ao%n)
!
   end subroutine destruct_frozen_CCT_hf
!
!
   module subroutine initialize_tdhf_quantities_hf(wf, n_states)
!!
!!    Initialize tdhf quantities
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      integer, intent(in) :: n_states
!
      wf%n_tdhf_vectors = n_states
!
      call wf%initialize_tdhf_files()
      call mem%alloc(wf%tdhf_excitation_energies, wf%n_tdhf_vectors)
!
   end subroutine initialize_tdhf_quantities_hf
!
!
   module subroutine destruct_tdhf_excitation_energies_hf(wf)
!!
!!    Destruct TDHF excitation energies
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      if (allocated(wf%tdhf_excitation_energies)) &
         call mem%dealloc(wf%tdhf_excitation_energies, wf%n_tdhf_vectors)
!
   end subroutine destruct_tdhf_excitation_energies_hf
!
!
end submodule initialize_destruct_hf
