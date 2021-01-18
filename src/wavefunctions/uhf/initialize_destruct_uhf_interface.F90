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
   module subroutine initialize_orbitals_uhf(wf)
!!
!!    Initialize orbitals
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the arrays associated with the orbital
!!    coefficients. In spin-unrestricted wavefunctions, this
!!    will include alpha and beta coefficients, though these
!!    are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine initialize_orbitals_uhf
!
!
   module subroutine initialize_density_uhf(wf)
!!
!!    Initialize density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the AO density (or densities).
!!    In spin-unrestricted wavefunctions, this alpha and beta densities,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine initialize_density_uhf
!
!
   module subroutine initialize_fock_uhf(wf)
!!
!!    Initialize Fock
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the AO Fock matrix (or matrices).
!!    In spin-unrestricted wavefunctions, this alpha and beta matrices,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine initialize_fock_uhf
!
!
   module subroutine destruct_fock_uhf(wf)
!!
!!    Destruct Fock
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Destructs the AO Fock matrix (or matrices).
!!    In spin-unrestricted wavefunctions, this alpha and beta matrices,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine destruct_fock_uhf
!
!
   module subroutine initialize_ao_density_a_uhf(wf)
!!
!!    Initialize AO density alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine initialize_ao_density_a_uhf
!
!
   module subroutine destruct_ao_density_uhf(wf)
!!
!!    Destruct AO density
!!    Written by Linda Goletto, Oct 2019
!!    Based on the work of Eirik F. Kjønstad, Sep 2018
!!
!!    Destructs the AO density matrix (or matrices).
!!    In spin-unrestricted wavefunctions, this alpha and beta matrices,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine destruct_ao_density_uhf
!
!
   module subroutine destruct_ao_density_a_uhf(wf)
!!
!!    Destruct AO density alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine destruct_ao_density_a_uhf
!
!
   module subroutine initialize_ao_density_b_uhf(wf)
!!
!!    Initialize AO density beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine initialize_ao_density_b_uhf
!
!
   module subroutine destruct_ao_density_b_uhf(wf)
!!
!!    Destruct AO density beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine destruct_ao_density_b_uhf
!
!
   module subroutine initialize_ao_fock_a_uhf(wf)
!!
!!    Initialize AO Fock alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine initialize_ao_fock_a_uhf
!
!
   module subroutine destruct_ao_fock_a_uhf(wf)
!!
!!    Destruct AO Fock alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine destruct_ao_fock_a_uhf
!
!
   module subroutine initialize_ao_fock_b_uhf(wf)
!!
!!    Initialize AO Fock beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine initialize_ao_fock_b_uhf
!
!
   module subroutine destruct_ao_fock_b_uhf(wf)
!!
!!    Destruct AO Fock beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine destruct_ao_fock_b_uhf
!
!
   module subroutine destruct_orbital_coefficients_uhf(wf)
!!
!!    Destruct orbital coefficients
!!    Written by Linda Goletto, Oct 2019
!!    Based on the work of Eirik F. Kjønstad, Sep 2018
!!
!!    Destructs the orbital coefficients matrix (or matrices).
!!    In spin-unrestricted wavefunctions, this alpha and beta matrices,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine destruct_orbital_coefficients_uhf
!
!
   module subroutine initialize_orbital_coefficients_a_uhf(wf)
!!
!!    Initialize orbital coefficients alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine initialize_orbital_coefficients_a_uhf
!
!
   module subroutine destruct_orbital_coefficients_a_uhf(wf)
!!
!!    Destruct orbital coefficients alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine destruct_orbital_coefficients_a_uhf
!
!
   module subroutine initialize_orbital_coefficients_b_uhf(wf)
!!
!!    Initialize orbital coefficients beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine initialize_orbital_coefficients_b_uhf
!
!
   module subroutine destruct_orbital_coefficients_b_uhf(wf)
!!
!!    Destruct orbital coefficients beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine destruct_orbital_coefficients_b_uhf
!
!
   module subroutine initialize_orbital_energies_a_uhf(wf)
!!
!!    Initialize orbital energies alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine initialize_orbital_energies_a_uhf
!
!
   module subroutine destruct_orbital_energies_uhf(wf)
!!
!!    Destruct orbital energies
!!    Written by Linda Goletto, Oct 2019
!!    Based on the work of Eirik F. Kjønstad, Sep 2018
!!
!!    Destructs the orbital energies array (or arrays).
!!    In spin-unrestricted wavefunctions, this alpha and beta arrays,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine destruct_orbital_energies_uhf
!
!
   module subroutine destruct_orbital_energies_a_uhf(wf)
!!
!!    Initialize orbital energies alpha
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine destruct_orbital_energies_a_uhf
!
!
   module subroutine initialize_orbital_energies_b_uhf(wf)
!!
!!    Initialize orbital energies beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine initialize_orbital_energies_b_uhf
!
!
   module subroutine destruct_orbital_energies_b_uhf(wf)
!!
!!    Destruct orbital energies beta
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine destruct_orbital_energies_b_uhf
