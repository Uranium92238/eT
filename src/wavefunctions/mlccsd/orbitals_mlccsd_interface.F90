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
!!    Based on which orbital types are requested on 
!!    input, we construct the new orbital basis 
!!    and set:
!!
!!       - n_ccsd_o, n_ccsd_v, n_cc2_o, n_cc2_v, n_ccs_o, n_ccs_v
!!
!!       - first_ccsd_o, last_ccsd_o, first_ccsd_v, last_ccsd_v
!!
!!       - first_cc2_o, last_cc2_o, first_cc2_v, last_cc2_v
!!
!!       - first_ccs_o, last_ccs_o, first_ccs_v, last_ccs_v
!!
!!    NOTE: This routine should always be 
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
!!    constructs the transformation matrices (occupied and virtual) 
!!    of that transformation.
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
   module subroutine construct_cholesky_orbitals_mlccsd(wf, occupied_only)
!!
!!    Construct Cholesky orbitals
!!    Written by Sarai D. Folkestad, Feb 2019
!!
!!    Constructs Cholesky orbitals by
!!    decomposing the HF AO density and the
!!    HF virtual AO density
!!
!!    See A. M. J. Sánchez de Merás, H. Koch, 
!!    I. G. Cuesta, and L. Boman (J. Chem. Phys. 132, 204105 (2010))
!!    for more information on active space generation
!!    using Cholesky decomposition
!!
!!    'occupied_only' : Optional argument that is used to 
!!                      determine if we construct occuped 
!!                      Cholesky orbitals only, or if we also 
!!                      construct the virtual Cholesky orbitals
!!                      Default: .false.
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
      logical, intent(in), optional :: occupied_only
!
   end subroutine construct_cholesky_orbitals_mlccsd
!
!
   module subroutine construct_paos_mlccsd(wf)
!!
!!    Construct PAOs
!!    Written by Sarai D. Folkestad, Feb 2019
!!
!!    Construct projected atomic orbitals for
!!    virtual orbitals.
!!
!!    1. Construct PAOs on active atoms
!!
!!    2. Orthonormalize these active virtual orbitals 
!!
!!    3. PAOs for the remaining system by projecting out
!!       both occupied and active virtuals out of all  AOs
!!
!!    4. Orthonormalize these inactive virtual orbitals 
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
   end subroutine construct_paos_mlccsd
!
!
   module subroutine construct_cc2_cnto_transformation_matrices_mlccsd(wf, T_o, T_v)
!!
!!    Construct CC2 CNTO transformation matrices
!!    Written by Sarai D. Folkestad, May 2019
!!
!!       - Run CCS calculation
!!
!!       - Construct M and N for CNTOs 
!!
!!       - Diagonalize M and N
!!
!!       - Write transformation matrices to file
!!
      implicit none
!
      class(mlccsd) :: wf
      real(dp), dimension(wf%n_o, wf%n_o), intent(out) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(out) :: T_v
!
   end subroutine construct_cc2_cnto_transformation_matrices_mlccsd
!
!
   module subroutine cc2_calculation_for_cntos_mlccsd(wf, transformation, n_cnto_states, &
                                                            R_ai, R_aibj, cnto_states)
!!
!!    CC2 calculation for CNTOs
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!
!!    Performs CC2 calculation for CNTO
!!    construction.
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
      character(len=200), intent(in) :: transformation
      integer, intent(in) :: n_cnto_states
      real(dp), dimension(wf%n_v, wf%n_o, n_cnto_states), intent(out) :: R_ai
      real(dp), dimension(wf%n_t1*(wf%n_t1+1)/2, n_cnto_states), intent(out) :: R_aibj
      integer, dimension(n_cnto_states), intent(in) :: cnto_states
!
   end subroutine cc2_calculation_for_cntos_mlccsd
