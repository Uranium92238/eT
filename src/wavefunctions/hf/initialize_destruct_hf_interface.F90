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
   module subroutine initialize_ao_density_hf(wf)
!!
!!    Initialize AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine initialize_ao_density_hf
!
!
   module subroutine initialize_ao_overlap_hf(wf)
!!
!!    Initialize AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine initialize_ao_overlap_hf
!
!
   module subroutine initialize_pivot_matrix_ao_overlap_hf(wf)
!!
!!    Initialize pivot matrix AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine initialize_pivot_matrix_ao_overlap_hf
!
!
   module subroutine initialize_cholesky_ao_overlap_hf(wf)
!!
!!    Initialize cholesky vectors AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine initialize_cholesky_ao_overlap_hf
!
!
   module subroutine initialize_shp_eri_schwarz_hf(wf)
!!
!!    Initialize shell pair eri schwarz
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine initialize_shp_eri_schwarz_hf
!
!
   module subroutine destruct_shp_eri_schwarz_hf(wf)
!!
!!    Destruct shell pair eri schwarz
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine destruct_shp_eri_schwarz_hf
!
!
   module subroutine initialize_shp_eri_schwarz_list_hf(wf)
!!
!!    Initialize shell pair eri schwarz list
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine initialize_shp_eri_schwarz_list_hf
!
!
   module subroutine destruct_shp_eri_schwarz_list_hf(wf)
!!
!!    Destruct shell pair eri schwarz list
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine destruct_shp_eri_schwarz_list_hf
!
!
   module subroutine destruct_ao_overlap_hf(wf)
!!
!!    Destruct AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine destruct_ao_overlap_hf
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
   end subroutine destruct_ao_density_hf
!
!
   module subroutine destruct_pivot_matrix_ao_overlap_hf(wf)
!!
!!    Destruct pivot matrix AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine destruct_pivot_matrix_ao_overlap_hf
!
!
   module subroutine destruct_cholesky_ao_overlap_hf(wf)
!!
!!    Initialize cholesky vectors AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine destruct_cholesky_ao_overlap_hf
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
   end subroutine destruct_fock_hf
!
!
   module subroutine initialize_ao_h_hf(wf)
!!
!!    Initialize AO h
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine initialize_ao_h_hf
!
!
   module subroutine destruct_ao_h_hf(wf)
!!
!!    Destruct AO h
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine destruct_ao_h_hf
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
   end subroutine destruct_orbital_coefficients_fc_hf
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
   end subroutine destruct_frozen_CCT_hf
