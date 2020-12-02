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
   module subroutine orbital_partitioning_mlcc2(wf)
!!
!!    Orbital partitioning 
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    This routine drives the orbital 
!!    partitioning for MLCC2 
!!
!!    Based on which orbitals are requested on 
!!    input, we construct this new orbital basis 
!!    and sets, n_cc2_o, n_cc2_v, n_ccs_o and n_ccs_v,
!!    first_cc2_o, last_cc2_o, first_cc2_v, last_cc2_v
!!    first_ccs_o, last_ccs_o, first_ccs_v, last_ccs_v
!!
!!    NOTE: This routine shold always be 
!!    followed by a routine which block diagonalizes
!!    the Fock matrix!
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
   end subroutine orbital_partitioning_mlcc2
!
!
   module subroutine construct_cholesky_orbitals_mlcc2(wf, occupied_only)
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
!!                      default: .false.
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
      logical, intent(in), optional :: occupied_only
!
   end subroutine construct_cholesky_orbitals_mlcc2
!
!
   module subroutine construct_block_diagonal_fock_orbitals_mlcc2(wf, n_levels, n_occupied_list, &
                                          n_virtual_list, orbital_coefficients, orbital_energies) 
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
      implicit none
!
      class(mlcc2), intent(inout) :: wf
      integer, intent(in) :: n_levels
      integer, dimension(n_levels), intent(in) :: n_occupied_list, n_virtual_list
      real(dp), dimension(wf%n_ao, wf%n_mo), intent(inout) :: orbital_coefficients
      real(dp), dimension(wf%n_mo), intent(inout) :: orbital_energies
!
   end subroutine construct_block_diagonal_fock_orbitals_mlcc2
!
!
   module subroutine construct_M_and_N_cnto_mlcc2(wf, R_ai, R_aibj, M, N, set_to_zero)
!!
!!    Construct M and N
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Constructs the M and N matrices,
!!
!!       M_ij += ( sum_a R_ai R_aj + 1/2 sum_abl(1 + δ_ai,bl δ_i,j) R_aibl R_ajbl )
!!       N_ab += ( sum_i R_ai R_bi + 1/2 sum_cij(1 + δ_ai,cj δ_a,b) R_aicj R_bicj )
!!
!!    Used to construct CNTOs.
!!
!!    set_to_zero determines if M and N are set to zero initially. This makes it possible for 
!!    the routine to be used to add to M and N if we are using more than one excitation vector 
!!    to generate the CNTOs.
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)  :: R_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                  :: R_ai
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: M
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: N
      logical, intent(in) :: set_to_zero
!
   end subroutine construct_M_and_N_cnto_mlcc2
!
!
   module subroutine construct_cntos_mlcc2(wf, T_o, T_v)
!!
!!    Construct correlated natural transition orbitals
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Transform orbital coefficients to CNTO
!!    basis.
!!
      implicit none
!
      class(mlcc2) :: wf
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: T_v
!
   end subroutine construct_cntos_mlcc2
!
!
   module subroutine read_cnto_transformation_matrices_mlcc2(wf, T_o, T_v)
!!
!!    Read CNTO transformation matrices
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Read CNTO transformation matrices.
!!
      implicit none
!
      class(mlcc2) :: wf
      real(dp), dimension(wf%n_o, wf%n_o), intent(out) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(out) :: T_v
!
   end subroutine read_cnto_transformation_matrices_mlcc2
!
!
   module subroutine write_cnto_transformation_matrices_mlcc2(wf, T_o, T_v)
!!
!!    Write CNTO transformation matrices
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Write CNTO transformation matrices.
!!    Used to ensure restart
!!
      implicit none
!
      class(mlcc2) :: wf
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: T_v
!
   end subroutine write_cnto_transformation_matrices_mlcc2
!
!
   module subroutine construct_ccs_cnto_transformation_matrices_mlcc2(wf, T_o, T_v)
!!
!!    Construct CCS CNTO transformation matrices
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
      class(mlcc2) :: wf
      real(dp), dimension(wf%n_o, wf%n_o), intent(out) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(out) :: T_v
!
   end subroutine construct_ccs_cnto_transformation_matrices_mlcc2
!
!
   module subroutine construct_M_nto_mlcc2(wf, R_ai, M, set_to_zero)
!!
!!    Construct M and N
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Constructs the M and N matrices,
!!
!!       M_ij +=  sum_a R_ai R_aj
!!
!!    Used to construct occupied NTOs.
!!
!!    set_to_zero determines if contributions are added to, or overwrites M and N. 
!!
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: M
      logical, intent(in) :: set_to_zero
!
   end subroutine construct_M_nto_mlcc2
!
!
   module subroutine read_nto_transformation_matrix_mlcc2(wf, T_o)
!!
!!    Read NTO transformation matrices
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Read NTO transformation matrices.
!!
      implicit none
!
      class(mlcc2) :: wf
      real(dp), dimension(wf%n_o, wf%n_o), intent(out) :: T_o
!
   end subroutine read_nto_transformation_matrix_mlcc2
!
!
   module subroutine write_nto_transformation_matrix_mlcc2(wf, T_o)
!!
!!    Write NTO transformation matrices
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Write NTO transformation matrices.
!!    Used to ensure restart
!!
      implicit none
!
      class(mlcc2) :: wf
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: T_o
!
   end subroutine write_nto_transformation_matrix_mlcc2
!
!
   module subroutine construct_ccs_nto_transformation_matrix_mlcc2(wf, T_o)
!!
!!    Construct CCS CNTO transformation matrices
!!    Written by Sarai D. Folkestad, May 2019
!!
!!       - Run CCS calculation
!!
!!       - Construct M NTOs 
!!
!!       - Diagonalize M
!!
!!       - Write transformation matrix to file
!!
      implicit none
!
      class(mlcc2) :: wf
      real(dp), dimension(wf%n_o, wf%n_o), intent(out) :: T_o
!
   end subroutine construct_ccs_nto_transformation_matrix_mlcc2
!
!
   module subroutine construct_mixed_nto_canonical_orbitals_mlcc2(wf, T_o)
!!
!!    Construct mixed NTO and canonical orbitals
!!    Written by Sarai D. Folekstad, Jun 2019
!!
!!    Construct occupiued NTOs, leave canonical virtuals
!!
      implicit none 
!
      class(mlcc2) :: wf
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: T_o
!
   end subroutine construct_mixed_nto_canonical_orbitals_mlcc2
!
!
   module subroutine construct_paos_mlcc2(wf)
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
      class(mlcc2), intent(inout) :: wf
!
   end subroutine construct_paos_mlcc2
!
!
   module subroutine diagonalize_M_and_N_mlcc2(wf, T_o, T_v)
!!
!!    Diagonalize M and N
!!    Written by Sarai D. Folkestad, Aug. 2019
!!
!!    Diagonalizes the M and N CNTO matrices.
!!
      implicit none
!
      class(mlcc2) :: wf
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: T_v
!
   end subroutine diagonalize_M_and_N_mlcc2
!
!
   module subroutine ccs_calculation_for_cntos_mlcc2(wf, transformation, n_cnto_states, &
                                                         R_ai, cnto_states, omega_ccs)
!!
!!    CCS calculation for CNTOs
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Performs CCS calculation for CNTO (and NTO)
!!    construction.
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
      character(len=200), intent(in) :: transformation
      integer, intent(in) :: n_cnto_states
      real(dp), dimension(wf%n_v, wf%n_o, n_cnto_states), intent(out) :: R_ai
      integer, dimension(n_cnto_states), intent(in) :: cnto_states
      real(dp), dimension(n_cnto_states), intent(out), optional :: omega_ccs
!
   end subroutine ccs_calculation_for_cntos_mlcc2
!
!
   module subroutine check_orthonormality_of_MOs_mlcc2(wf)
!!
!!    Check orthonormality of MOs
!!    Write Sarai D. Folkestad, Sep 2019
!!
!!    Checks that 
!!
!!       C^T S C = I
!!
!!    to ensure that we have orthonormal MOs.
!!
      implicit none
!
      class(mlcc2) :: wf
!
   end subroutine check_orthonormality_of_MOs_mlcc2
!
!
   module subroutine add_doubles_M_and_N_cnto_mlcc2(wf, M, N, doubles_file)
!!
!!    Add doubles M and N CNTO 
!!    Written by Sarai D. Folkestad, Feb 2020
!!
!!    The CNTO matrices are defined as
!!
!!       M_ij = ( sum_a R_ai R_aj + 1/2 sum_abl(1 + δ_ai,bl δ_i,j) R_aibl R_ajbl )
!!       N_ab = ( sum_i R_ai R_bi + 1/2 sum_cij(1 + δ_ai,cj δ_a,b) R_aicj R_bicj )
!!
!!    In this routine the doubles part 
!!
!!       M_ij += ( 1/2 sum_abl(1 + δ_ai,bl δ_i,j) R_aibl R_ajbl )
!!       N_ab += ( 1/2 sum_cij(1 + δ_ai,cj δ_a,b) R_aicj R_bicj )
!!
!!    is added for the R_aibj on doubles_file
!!
      implicit none
!
      class(mlcc2) :: wf
      type(direct_stream_file) :: doubles_file
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: M
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: N
!
   end subroutine add_doubles_M_and_N_cnto_mlcc2
!
!
   module subroutine construct_M_and_N_singles_cnto_mlcc2(wf, R_ai, M, N, set_to_zero)
!!
!!    Construct M and N
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Constructs the M and N matrices,
!!
!!       M_ij = ( sum_a R_ai R_aj )
!!       N_ab = ( sum_i R_ai R_bi )
!!
!!    Used to construct CNTOs.
!!
!!    set_to_zero determines if contributions are added to, or overwrites M and N. 
!!
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: M
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: N
      logical, intent(in) :: set_to_zero
!
   end subroutine construct_M_and_N_singles_cnto_mlcc2
!
!
   module subroutine construct_semicanonical_mlcc_orbitals_mlcc2(wf)
!!
!!    Construct semicanonical orbitals
!!    Written by Sarai D. Folkestad
!!
!!    Wrapper to construct orbitals that block diagonalizes the fock matrix
!!    for the different orbitals
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
   end subroutine construct_semicanonical_mlcc_orbitals_mlcc2
