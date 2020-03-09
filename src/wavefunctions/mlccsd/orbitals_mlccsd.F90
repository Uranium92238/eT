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
submodule (mlccsd_class) orbitals_mlccsd
!
!!
!!    MLCCSD orbitals
!!
!!
!!    This submodule contains routines that handle orbital
!!    transformation and orbital partitioning for MLCCSD.
!!
!!    Orbital construction in MLCC calculations consists of 3
!!    steps:
!!
!!    1. Construct the orbital basis specified on input
!!
!!    2. If the number of active/inactive orbitals not determined
!!       on input, it is determined during orbital construction.
!!       Orbital coefficient matrix is ordered according to 
!!       active and inactive:
!!       
!!          C = [active occupied, inactive occupied, active virtual, inactive virtual]
!!
!!    3. The Fock matrix is constructed, and the F_oo and F_vv blocks are block 
!!       diagonalized, such that the active-active and inactive-inactive
!!       blocks are diagonal. The corresponding basis is the basis used in the 
!!       MLCC calculation.
!!
!!    Note that step 1 and 2 happens at the same time, in the routines 
!!    which constructs the orbitals
!!
!
   implicit none
!
contains
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
      use eri_cd_class, only : eri_cd
!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      type(timings) :: timer
!
      real(dp), dimension(:,:), allocatable :: T_o, T_v ! CNTO transformation matrices
!
      timer = timings('MLCC orbital construction', pl='minimal')
      call timer%turn_on()
!
      if (trim(wf%ccsd_orbital_type) == 'cnto-approx') then
!
!        Set orbital space sizes
!
         if (wf%do_cc2 .and. .not. wf%do_ccs) then
!
            wf%n_cc2_o = wf%n_o - wf%n_ccsd_o
            wf%n_cc2_v = wf%n_v - wf%n_ccsd_v
!  
         elseif (.not. wf%do_cc2) then
!
            wf%n_cc2_o = 0
            wf%n_cc2_v = 0
!
         endif
!
         if (wf%do_ccs) then
!
            wf%n_ccs_o = wf%n_o - wf%n_cc2_o - wf%n_ccsd_o
            wf%n_ccs_v = wf%n_v - wf%n_cc2_v - wf%n_ccsd_v
!
         else
!
            wf%n_ccs_o = 0
            wf%n_ccs_v = 0
!
         endif
!
         call mem%alloc(T_o, wf%n_o, wf%n_o)
         call mem%alloc(T_v, wf%n_v, wf%n_v)
!
         if (wf%cnto_restart) then
!
            call wf%read_cnto_transformation_matrices(T_o, T_v)
!
         else
!
            if (trim(wf%ccsd_orbital_type) == 'cnto-approx') then
!
               call wf%construct_ccs_cnto_transformation_matrices(T_o, T_v)
!
            endif
!
         endif
!
         call wf%construct_cntos(T_o, T_v)
!
         call mem%dealloc(T_o, wf%n_o, wf%n_o)
         call mem%dealloc(T_v, wf%n_v, wf%n_v)
!
      else
!
         call output%error_msg('Could not recognize the specified orbital type')
!
      endif
!
!     Set orbital partitioning specifications
!
      wf%first_cc2_o = wf%n_ccsd_o + 1
      wf%first_cc2_v = wf%n_ccsd_v + 1
!
      wf%last_cc2_o = wf%first_cc2_o + wf%n_cc2_o - 1
      wf%last_cc2_v = wf%first_cc2_v + wf%n_cc2_v - 1
!
      wf%first_ccs_o = wf%last_cc2_o + 1
      wf%first_ccs_v = wf%last_cc2_v + 1
!
      wf%last_ccs_o = wf%n_o
      wf%last_ccs_v = wf%n_v
!
      call timer%turn_off()
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
      real(dp), dimension(:,:), allocatable :: S
      real(dp), dimension(:,:), allocatable :: X
!
      integer :: n_doubles_v, n_doubles_o
!
      n_doubles_o = wf%n_ccsd_o + wf%n_cc2_o
      n_doubles_v = wf%n_ccsd_v + wf%n_cc2_v
!
      call mem%alloc(S, wf%n_ao, wf%n_ao)
      call wf%get_ao_s_wx(S)
!
      call mem%alloc(X, n_doubles_o, wf%n_ao)
!
      call dgemm('T', 'N',                      &
                  n_doubles_o,                  &
                  wf%n_ao,                      &
                  wf%n_ao,                      &
                  one,                          &
                  wf%orbital_coefficients_cc2,  &
                  wf%n_ao,                      &
                  S,                            &
                  wf%n_ao,                      &
                  zero,                         &
                  X,                            &
                  n_doubles_o)
!
      call dgemm('N', 'N',                &
                  n_doubles_o,            &
                  n_doubles_o,            &
                  wf%n_ao,                &
                  one,                    &
                  X,                      &
                  n_doubles_o,            &
                  wf%orbital_coefficients,&
                  wf%n_ao,                &
                  zero,                   &
                  wf%O_o,                 &
                  n_doubles_o)
!
      call mem%dealloc(X, n_doubles_o, wf%n_ao)
!
      call mem%alloc(X, n_doubles_v, wf%n_ao)
!
      call dgemm('T', 'N',                                      &
                  n_doubles_v,                                  &
                  wf%n_ao,                                      &
                  wf%n_ao,                                      &
                  one,                                          &
                  wf%orbital_coefficients_cc2(1, wf%n_o + 1),   &
                  wf%n_ao,                                      &
                  S,                                            &
                  wf%n_ao,                                      &
                  zero,                                         &
                  X,                                            &
                  n_doubles_v)
!
      call mem%dealloc(S, wf%n_ao, wf%n_ao)
!
      call dgemm('N', 'N',                               &
                  n_doubles_v,                           &
                  n_doubles_v,                           &
                  wf%n_ao,                               &
                  one,                                   &
                  X,                                     &
                  n_doubles_v,                           &
                  wf%orbital_coefficients(1, wf%n_o + 1),&
                  wf%n_ao,                               &
                  zero,                                  &
                  wf%O_v,                                &
                  n_doubles_v)
!
      call mem%dealloc(X, n_doubles_v, wf%n_ao)
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
!
      use array_utilities, only : block_diagonalization
! 
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      real(dp), dimension(:,:), allocatable :: F_oo, F_vv, C_copy
!
      call mem%alloc(F_oo, wf%n_o, wf%n_o)
      call mem%alloc(F_vv, wf%n_v, wf%n_v)
!
      call dcopy(wf%n_o**2, wf%fock_ij, 1, F_oo, 1)
      call dcopy(wf%n_v**2, wf%fock_ab, 1, F_vv, 1)
!
      if (wf%do_cc2) then
!
         call dcopy(wf%n_mo*wf%n_ao, wf%orbital_coefficients, 1, wf%orbital_coefficients_cc2, 1)
         call dcopy(wf%n_mo, wf%orbital_energies, 1, wf%orbital_energies_cc2, 1)
!
      endif
!
!     Block diagonal occupied-occupied Fock
!
      call block_diagonalization(F_oo, wf%n_o, 3, &
                                 [integer::wf%n_ccsd_o, wf%n_cc2_o, wf%n_ccs_o], &
                                 wf%orbital_energies(1:wf%n_o))
!
      call block_diagonalization(F_vv, wf%n_v, 3, &
                                 [integer::wf%n_ccsd_v, wf%n_cc2_v, wf%n_ccs_v], &
                                 wf%orbital_energies(wf%n_o + 1 : wf%n_mo))
!
!     Transform blocks
!
      call mem%alloc(C_copy, wf%n_ao, wf%n_mo)
      call dcopy(wf%n_mo*wf%n_ao, wf%orbital_coefficients, 1, C_copy, 1)
      call zero_array(wf%orbital_coefficients, wf%n_mo*wf%n_ao)
!
!     1. CCSD occcupied
!
      if (wf%n_ccsd_o .gt. 0) then
!
         call dgemm('N', 'N',                 &
                     wf%n_ao,                 &
                     wf%n_ccsd_o,             &
                     wf%n_ccsd_o,             &
                     one,                     &
                     C_copy,                  &
                     wf%n_ao,                 &
                     F_oo,                    &
                     wf%n_o,                  &
                     one,                     &
                     wf%orbital_coefficients, &
                     wf%n_ao)
!
      endif
!
!     2. CC2 occcupied
!
      if (wf%n_cc2_o .gt. 0) then
!
         call dgemm('N', 'N',                                     &
                     wf%n_ao,                                     &
                     wf%n_cc2_o,                                  &
                     wf%n_cc2_o,                                  &
                     one,                                         &
                     C_copy(1, wf%n_ccsd_o + 1),                  &
                     wf%n_ao,                                     &
                     F_oo(wf%n_ccsd_o + 1, wf%n_ccsd_o + 1),      &
                     wf%n_o,                                      &
                     one,                                         &
                     wf%orbital_coefficients(1, wf%n_ccsd_o + 1), &
                     wf%n_ao)
      endif
!
!     3. CCS occupied
!
      if ((wf%n_ccsd_o + wf%n_cc2_o) .lt. wf%n_o) then
!
         call dgemm('N', 'N',                                                          &
                     wf%n_ao,                                                          &
                     wf%n_ccs_o,                                                       &
                     wf%n_ccs_o,                                                       &
                     one,                                                              &
                     C_copy(1, wf%n_ccsd_o + wf%n_cc2_o + 1),                          &
                     wf%n_ao,                                                          &
                     F_oo(wf%n_ccsd_o + wf%n_cc2_o + 1, wf%n_ccsd_o + wf%n_cc2_o + 1), &
                     wf%n_o,                                                           &
                     one,                                                              &
                     wf%orbital_coefficients(1, wf%n_ccsd_o + wf%n_cc2_o + 1),         &
                     wf%n_ao)
      endif
!
!     4. CCSD virtual
!
      if (wf%n_ccsd_v .gt. 0) then
!
         call dgemm('N', 'N',                               &
                     wf%n_ao,                               &
                     wf%n_ccsd_v,                           &
                     wf%n_ccsd_v,                           &
                     one,                                   &
                     C_copy(1, wf%n_o + 1),                 &
                     wf%n_ao,                               &
                     F_vv,                                  &
                     wf%n_v,                                &
                     one,                                   &
                     wf%orbital_coefficients(1, wf%n_o + 1),&
                     wf%n_ao)
!
      endif
!
!     5. CC2 virtual
!
      if (wf%n_cc2_v .gt. 0) then
!
         call dgemm('N', 'N',                                              &
                     wf%n_ao,                                              &
                     wf%n_cc2_v,                                           &
                     wf%n_cc2_v,                                           &
                     one,                                                  &
                     C_copy(1, wf%n_o + wf%n_ccsd_v + 1),                  &
                     wf%n_ao,                                              &
                     F_vv(wf%n_ccsd_v + 1, wf%n_ccsd_v + 1),               &
                     wf%n_v,                                               &
                     one,                                                  &
                     wf%orbital_coefficients(1, wf%n_o + wf%n_ccsd_v + 1), &
                     wf%n_ao)
      endif
!
!     6. CCS virtual
!
      if ((wf%n_ccsd_v + wf%n_cc2_v) .lt. wf%n_v) then
!
         call dgemm('N', 'N',                                                          &
                     wf%n_ao,                                                          &
                     wf%n_ccs_v,                                                       &
                     wf%n_ccs_v,                                                       &
                     one,                                                              &
                     C_copy(1, wf%n_o + wf%n_ccsd_v + wf%n_cc2_v + 1),                 &
                     wf%n_ao,                                                          &
                     F_vv(wf%n_ccsd_v + wf%n_cc2_v + 1, wf%n_ccsd_v + wf%n_cc2_v + 1), &
                     wf%n_v,                                                           &
                     one,                                                              &
                     wf%orbital_coefficients(1, wf%n_o + wf%n_ccsd_v + wf%n_cc2_v + 1),&
                     wf%n_ao)
      endif
!
      if (wf%do_cc2) then
!  
         call dcopy(wf%n_o**2, wf%fock_ij, 1, F_oo, 1)
         call dcopy(wf%n_v**2, wf%fock_ab, 1, F_vv, 1)
!
!        Block diagonal occupied-occupied Fock
!
         call block_diagonalization(F_oo, wf%n_o, 2, &
                                    [integer::wf%n_ccsd_o + wf%n_cc2_o, wf%n_ccs_o], &
                                    wf%orbital_energies_cc2(1:wf%n_o))
!
!        Block diagonal virtual-virtual Fock
!
         call block_diagonalization(F_vv, wf%n_v, 2, [&
                                    integer::wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccs_v], &
                                    wf%orbital_energies_cc2(wf%n_o+1:wf%n_mo))
!
!        Transform blocks
!
         call dcopy(wf%n_mo*wf%n_ao, wf%orbital_coefficients_cc2, 1, C_copy, 1)
         call zero_array(wf%orbital_coefficients_cc2, wf%n_mo*wf%n_ao)
!
!        1. Active occcupied
!
         if ((wf%n_cc2_o + wf%n_ccsd_o) .gt. 0) then
!
            call dgemm('N', 'N',                      &
                        wf%n_ao,                      &
                        wf%n_cc2_o + wf%n_ccsd_o,     &
                        wf%n_cc2_o + wf%n_ccsd_o,     &
                        one,                          &
                        C_copy,                       &
                        wf%n_ao,                      &
                        F_oo,                         &
                        wf%n_o,                       &
                        one,                          &
                        wf%orbital_coefficients_cc2,  &
                        wf%n_ao)
!
         endif
!
!        2. Inactive occupied
!
         if ((wf%n_cc2_o + wf%n_ccsd_o) .lt. wf%n_o) then
!
            call dgemm('N', 'N',                                                          &
                        wf%n_ao,                                                          &
                        wf%n_ccs_o,                                                       &
                        wf%n_ccs_o,                                                       &
                        one,                                                              &
                        C_copy(1, wf%n_ccsd_o + wf%n_cc2_o + 1),                          &
                        wf%n_ao,                                                          &
                        F_oo(wf%n_ccsd_o + wf%n_cc2_o + 1, wf%n_ccsd_o + wf%n_cc2_o + 1), &
                        wf%n_o,                                                           &
                        one,                                                              &
                        wf%orbital_coefficients_cc2(1, wf%n_ccsd_o + wf%n_cc2_o + 1),     &
                        wf%n_ao)
         endif
!
!        3. Active virtual
!
         if (wf%n_ccsd_v + wf%n_cc2_v .gt. 0) then
!
            call dgemm('N', 'N',                                     &
                        wf%n_ao,                                     &
                        wf%n_ccsd_v + wf%n_cc2_v,                    &
                        wf%n_ccsd_v + wf%n_cc2_v,                    &
                        one,                                         &
                        C_copy(1, wf%n_o + 1),                       &
                        wf%n_ao,                                     &
                        F_vv,                                        &
                        wf%n_v,                                      &
                        one,                                         &
                        wf%orbital_coefficients_cc2(1, wf%n_o + 1),  &
                        wf%n_ao)
!
         endif
!
!        4. Inactive virtual
!
         if ((wf%n_ccsd_v + wf%n_cc2_v) .lt. wf%n_v) then
!
            call dgemm('N', 'N',                                                                &
                        wf%n_ao,                                                                &
                        wf%n_ccs_v,                                                             &
                        wf%n_ccs_v,                                                             &
                        one,                                                                    &
                        C_copy(1, wf%n_o + wf%n_ccsd_v + wf%n_cc2_v + 1),                       &
                        wf%n_ao,                                                                &
                        F_vv(wf%n_ccsd_v + wf%n_cc2_v + 1, wf%n_ccsd_v + wf%n_cc2_v + 1),       &
                        wf%n_v,                                                                 &
                        one,                                                                    &
                        wf%orbital_coefficients_cc2(1, wf%n_o + wf%n_ccsd_v + wf%n_cc2_v + 1),  &
                        wf%n_ao)
!
         endif
!
      endif
!
      call mem%dealloc(C_copy, wf%n_ao, wf%n_mo)
!
      call mem%dealloc(F_oo, wf%n_o, wf%n_o)
      call mem%dealloc(F_vv, wf%n_v, wf%n_v)
!
   end subroutine construct_block_diagonal_fock_orbitals_mlccsd
!
!
end submodule orbitals_mlccsd
