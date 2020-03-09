!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
   module subroutine orbital_partitioning_mlccsd(wf)
!!
!!    Orbital partitioning 
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    This routine drives the orbital 
!!    partitioning for MLCCSD
!!
!!    Based on which orbitals are requested on 
!!    input, we construct this new orbital basis 
!!    and sets:
!!
!!       - n_ccsd_o, n_ccsd_v, n_cc2_o, n_cc2_v, n_ccs_o, n_ccs_v
!!
!!       - first_ccsd_o, last_ccsd_o, first_ccsd_v, last_ccsd_v
!!
!!       - first_cc2_o, last_cc2_o, first_cc2_v, last_cc2_v
!!
!!       - first_ccs_o, last_ccs_o, first_ccs_v, last_ccs_v
!!
!!    NOTE: This routine shold always be 
!!    followed by a routine which block diagonalizes
!!    the Fock matrix!
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
   end subroutine orbital_partitioning_mlccsd
!
!
   module subroutine construct_mlccsd_basis_transformation_matrix_mlccsd(wf)
!!
!!    Construct mlccsd basis transformtion matrix
!!    Written by Sarai D. Folkestad, 2017-2019
!!
!!    For an MLCCSD calculation with CC2, we must swap MO basis
!!    to detemine the CC2 double amplitudes. This routine
!!    constructs the transfromation matrices (occupied and virtual) 
!!    of that transfromation.
!!
!!    The orthogonal transformation matrix
!!    is defined as
!!
!!       O = C1^T S C2, 
!!
!!    where C1 is the CC2 basis, C2 is the CCSD basis,
!!    and S is the AO overlap matrix.
!!
!!    Since this transformation only will mix occupied 
!!    orbitals with occupied orbitals, and virtual orbitals 
!!    with virtual orbitals, we store only the two diagonal blocks
!!
!!       O_o (n_o x n_o)
!!
!!    and 
!!
!!       O_v (n_v x n_v).
!!
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
   end subroutine construct_mlccsd_basis_transformation_matrix_mlccsd
!
!
   module subroutine construct_block_diagonal_fock_orbitals_mlccsd(wf)
!!
!!    Construct block diagonal Fock MOs 
!!    Written by Sarai D. Folkestad, Feb 2019
!!
!!    This routine constructs the MOs which
!!    block-diagonalize the occupied-occupied
!!    and virutal-virtual blocks of the Fock matrix s.t.
!!    the active-active, and inactive-inactive blocks are 
!!    diagonal.
!!   
!!    Note that after this routine, the Fock matrix in wf 
!!    corresponds to the old basis but the MOs are updated.
!!
!!    In the case of do_cc2 = .true. we must also update
!!    the CC2 orbital basis, such that it corresponds to 
!!    a block diagonal fock matrix, where the diagonal blocks 
!!    corresponding to CC2 and CCSD orbitals are diagonal.
!!    
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
   end subroutine construct_block_diagonal_fock_orbitals_mlccsd
!
!
   module subroutine construct_block_diagonal_fock_mos_3_level_mlccsd(wf, n_total, n_l_1, &
                                                                     n_l_2, fock, MO_coeff, diagonal)
!!
!!    Construct block diagonal fock 3 levels
!!    Written by Sarai D. Folkestad, Feb 2019
!!
!!    Construct orbitals that block diagonalizes 
!!    Fock matrix block for three levels (inactive/active)
!!
!!    'n_total' : Total dimmension of Fock matrix block
!!
!!    'n_active' : Dimension of active block of Fock matrix block
!!
!!    'fock' : Fock matrix block to block diagonalize
!!
!!    'mo_coef' : MO coefficients which are updated to 
!!                the new basis which block diagonalizes Fock
!!                matrix block
!!
!!    'diagonal' : On exit the, diagonal elements of Fock matrix after block
!!                 diagonalization
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
      integer, intent(in) :: n_total, n_l_1, n_l_2 ! Total matrix dimension of block, and number of active orbitals
      real(dp), dimension(n_total, n_total), intent(inout) :: fock
      real(dp), dimension(wf%n_ao, n_total), intent(inout) :: mo_coeff
      real(dp), dimension(n_total), intent(out)            :: diagonal
!
   end subroutine construct_block_diagonal_fock_mos_3_level_mlccsd
