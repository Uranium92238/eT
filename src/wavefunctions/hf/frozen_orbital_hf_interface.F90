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
   module subroutine prepare_mos_hf(wf)
!!
!!    Prepare MOs
!!    Written by Ida-Marie Høyvik, Oct 2019
!!
!!    This routine prepares the MOs for coupled cluster
!!    in the cases where there is a reduction in the
!!    number of MOs in CC compared to HF
!!
!!    Examples of this is the frozen core
!!    approximation and if CC is only done
!!    for a localized region of a large molecule
!!    which has been treated at HF level of theory.
!!
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine prepare_mos_hf
!
!
   module subroutine remove_core_orbitals_hf(wf)
!!
!!    Remove core orbitals
!!    Written by Sarai D. Folkestad, Sep 2018
!!
!!    - Removes core orbitals from wf%orbital_coefficients
!!
!!       Removes 1s for C - Mg
!!       Removes 1s, 2s, 2p for Al - Zn
!!
!!    - The core orbitals are stored in wf%orbital_coefficients_fc
!!    - The number of frozen core orbitals is wf%n_frozen_core_orbitals 
!!       on exit
!!    - On exit wf%n_mo and wf%n_o are updated to not include the core orbitals
!!    
      implicit none
!
      class(hf) :: wf
!
   end subroutine remove_core_orbitals_hf
!
!
   module subroutine remove_frozen_hf_orbitals_hf(wf)
!!
!!    Remove frozen hf orbitals
!!    Written by Sarai D. Folkestad, Feb 2019
!!    Added and modified for HF by Ida-Marie Hoyvik, Oct 2019
!!
!!    - Removes frozen HF molecular orbitals from wf%orbital_coefficients
!!      and places them in wf%orbital_coefficients_frozen_hf
!!    - The number of frozen occupied HF orbitals is wf%n_frozen_hf_o
!!      on exit
!!    - On exit wf%n_mo, wf%n_o and wf%n_v are updated not to include 
!!      frozen hf orbitals
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
   end subroutine remove_frozen_hf_orbitals_hf
!
!
   module subroutine diagonalize_fock_frozen_hf_orbitals_hf(wf)
!!
!!    Diagonalize Fock frozen HF orbitals
!!    Written by Sarai D. Folkestad and Linda Goletto, Nov 2019
!!
!!    Does a diagonalization of the Fock matrix in the 
!!    MO basis where the frozen HF orbitals have been removed
!!
!!    Fock matrix is no longer diagonal, because determining
!!    the frozen HF orbitals entails mixing of occupied orbitals 
!!    and mixing of virtual orbitals, respectively.
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine diagonalize_fock_frozen_hf_orbitals_hf
!
!
   module subroutine prepare_frozen_fock_terms_hf(wf)
!!
!!    Prepare frozen Fock contributions
!!    Written by Sarai D. Folkestad, Oct 2019
!!
!!    This routine prepares the frozen Fock contributions
!!    to coupled cluster. This occurs e.g.,  in the cases where there
!!    is a reduction in the number of MOs in CC compared to HF
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine prepare_frozen_fock_terms_hf
!
!
   module subroutine construct_mo_fock_fc_term_hf(wf)
!!
!!    Calculate MO Fock frozen core contribution
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Constructs the frozen core contribution to
!!    the fock matrix
!!
!!       F_pq = (2 g_wxyz D^FC_yz - g_wyzx D^FC_xy) C_pw C_qx
!!
!!    in preparation of FC-CC
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine construct_mo_fock_fc_term_hf
!
!
   module subroutine construct_mo_fock_frozen_hf_term_hf(wf)
!!
!!    Construct MO fock frozen hf  contribution
!!    Written by Ida-Marie Høyvik, Oct 2019
!!
!!
!!    Constructs the frozen HF contribution to
!!    the fock matrix
!!
!!       F_pq = (2 g_wxyz D^F_yz - g_wyzx D^F_xy) C_pw C_qx
!!
!!    in preparation of CC in subspace.
!!
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine construct_mo_fock_frozen_hf_term_hf
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
   module function get_n_active_hf_atoms_hf(wf) result(n_active_hf_atoms)
!!
!!    Get number of active hf atoms
!!    Written by Sarai D. Folkestad and Linda Goletto, Dec 2019
!!
!!    Sets the number of active hf atoms in the system 
!!
!!    NOTE: modified in mlhf
!!
      implicit none 
!
      class(hf), intent(in) :: wf
      integer :: n_active_hf_atoms
!
   end function get_n_active_hf_atoms_hf
