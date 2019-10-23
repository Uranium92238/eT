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
!
!
   module subroutine orbital_partitioning_mlcc2(wf)
!!
!!    Orbital partitioning 
!!    Written by Sarai D. Folkestad, May 2019
!!
      use eri_cd_class
!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
   end subroutine orbital_partitioning_mlcc2
!
!
   module subroutine construct_cholesky_orbitals_mlcc2(wf, occupied_only)
!!
!!    Cholesky orbitals
!!    Written by Sarai D. Folkestad, Feb 2019
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      logical, intent(in), optional :: occupied_only
!
   end subroutine construct_cholesky_orbitals_mlcc2
!
!
   module subroutine construct_block_diagonal_fock_mos_2_level_mlcc2(wf, n_total, n_active, fock, MO_coeff, diagonal)
!!
!!    Construct Fock block diagonal 2 levels
!!    Written by Sarai D. Folkestad, Feb 2019
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      integer, intent(in) :: n_total, n_active ! Total matrix dimension of block, and number of active orbitals
!
      real(dp), dimension(n_total, n_total), intent(inout) :: fock
      real(dp), dimension(n_total, n_total), intent(inout) :: mo_coeff
      real(dp), dimension(n_total), intent(out) :: diagonal
!
   end subroutine construct_block_diagonal_fock_mos_2_level_mlcc2
!
!
   module subroutine construct_block_diagonal_fock_orbitals_mlcc2(wf)
!!
!!    Construct orbitals that block diagonalizes Fock
!!    Written by Sarai D. Folkestad, Feb 2019
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
   end subroutine construct_block_diagonal_fock_orbitals_mlcc2
!
!
   module subroutine construct_M_and_N_cnto_mlcc2(wf, R_ai, R_aibj, M, N, set_to_zero)
!!
!!    Construct M and N
!!    Written by Sarai D. Folkestad, May 2019
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: R_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                 :: R_ai
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: M
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: N
!
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
      implicit none
!
      class(mlcc2) :: wf
!
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
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(out) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(out) :: T_v
!
   end subroutine read_cnto_transformation_matrices_mlcc2
!
!
   module subroutine construct_ccs_cnto_transformation_matrices_mlcc2(wf, T_o, T_v)
!!
!!    Construct CCS CNTO transformation matrices
!!    Written by Sarai D. Folkestad, May 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(out) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(out) :: T_v
!
   end subroutine construct_ccs_cnto_transformation_matrices_mlcc2

!
   module subroutine construct_M_nto_mlcc2(wf, R_ai, M, set_to_zero)
!!
!!    Construct M and N
!!    Written by Sarai D. Folkestad, May 2019
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: M
!
      logical, intent(in) :: set_to_zero
!
   end subroutine construct_M_nto_mlcc2
!
!
   module subroutine construct_ccs_nto_transformation_matrix_mlcc2(wf, T_o)
!!
!!    Construct CCS CNTO transformation matrices
!!    Written by Sarai D. Folkestad, May 2019
!!
!
      implicit none
!
      class(mlcc2) :: wf
!
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
!
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
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
   end subroutine construct_paos_mlcc2
!
!
   module subroutine diagonalize_M_and_N_mlcc2(wf,T_o, T_v)
!!
!!    Diagonalize M and N
!!    Written by Sarai D. Folkestad, Aug. 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
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
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      character(len=200), intent(in) :: transformation
!
      integer, intent(in) :: n_cnto_states
!
      real(dp), dimension(wf%n_v, wf%n_o, n_cnto_states), intent(out) :: R_ai
!
      integer, dimension(n_cnto_states), intent(in) :: cnto_states
!
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
!! 
      implicit none

      class(mlcc2) :: wf

   end subroutine check_orthonormality_of_MOs_mlcc2