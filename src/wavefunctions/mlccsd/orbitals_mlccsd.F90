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
!     Sanity check:  CC2 orbitals should be the same as CCSD orbitals
!                    for 3 level calculations.
!
      if (wf%do_ccs .and. wf%do_cc2) then
!
         if (wf%cc2_orbital_type .ne. wf%ccsd_orbital_type)  then
!
               call output%error_msg('CCSD and CC2 orbital types must be the same!')
         endif   
!  
      endif
!
      if (trim(wf%ccsd_orbital_type) == 'cnto-approx' .or. &
         trim(wf%ccsd_orbital_type) == 'cnto' ) then
!
!        Set orbital space sizes
!
!        Note: If do_ccs and do_cc2, then n_cc2_o and n_cc2_v 
!              is set directly from input.<
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
            call output%printf('m', 'Requested restart for CNTOs, &
                                    &reading orbital transformation matrices')
!
            call wf%read_cnto_transformation_matrices(T_o, T_v)
!
         else
!
            if (trim(wf%ccsd_orbital_type) == 'cnto-approx') then
!
               call wf%construct_ccs_cnto_transformation_matrices(T_o, T_v)
!
            elseif (trim(wf%ccsd_orbital_type) == 'cnto') then
!
               call wf%construct_cc2_cnto_transformation_matrices(T_o, T_v)
!      
            endif
!
            call wf%write_cnto_transformation_matrices(T_o, T_v)
!
         endif
!
         call wf%construct_cntos(T_o, T_v)
!
         call mem%dealloc(T_o, wf%n_o, wf%n_o)
         call mem%dealloc(T_v, wf%n_v, wf%n_v)
!
      elseif (trim(wf%ccsd_orbital_type) == 'cholesky') then
!
         call wf%construct_cholesky_orbitals()
!
      elseif (trim(wf%ccsd_orbital_type) == 'cholesky-pao') then
!
         call wf%construct_cholesky_orbitals(occupied_only=.true.)
!
         wf%n_cc2_v = 0
         wf%n_ccs_v = 0
         wf%n_ccsd_v = 0
!
         call wf%construct_paos()
!
      elseif (trim(wf%ccsd_orbital_type) == 'nto-canonical') then
!
         call mem%alloc(T_o, wf%n_o, wf%n_o)
!
         if (wf%nto_restart) then
!
            call output%printf('m', 'Requested restart for NTOs, &
                                    &reading orbital transformation matrix')
!
            call wf%read_nto_transformation_matrix(T_o)
!
         else

            call wf%construct_ccs_nto_transformation_matrix(T_o)
!
         endif
!
         call wf%construct_mixed_nto_canonical_orbitals(T_o)
!
         if ( .not. wf%do_cc2) then
!
            wf%n_cc2_o = 0
            wf%n_cc2_v = 0
!
         elseif (.not. wf%do_ccs) then
!
            wf%n_cc2_o = wf%n_o - wf%n_ccsd_o
            wf%n_cc2_v = wf%n_v - wf%n_ccsd_v
!
         endif
!
         wf%n_ccs_o = wf%n_o - wf%n_ccsd_o - wf%n_cc2_o
         wf%n_ccs_v = wf%n_v - wf%n_ccsd_v - wf%n_cc2_v
!
         call mem%dealloc(T_o, wf%n_o, wf%n_o)
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
      real(dp), dimension(:,:), allocatable :: S
      real(dp), dimension(:,:), allocatable :: X
!
      integer :: n_doubles_v, n_doubles_o
!
      n_doubles_o = wf%n_ccsd_o + wf%n_cc2_o
      n_doubles_v = wf%n_ccsd_v + wf%n_cc2_v
!
      call mem%alloc(S, wf%ao%n, wf%ao%n)
      call wf%ao%get_oei('overlap', S)
!
      call mem%alloc(X, n_doubles_o, wf%ao%n)
!
      call dgemm('T', 'N',                      &
                  n_doubles_o,                  &
                  wf%ao%n,                      &
                  wf%ao%n,                      &
                  one,                          &
                  wf%orbital_coefficients_cc2,  &
                  wf%ao%n,                      &
                  S,                            &
                  wf%ao%n,                      &
                  zero,                         &
                  X,                            &
                  n_doubles_o)
!
      call dgemm('N', 'N',                &
                  n_doubles_o,            &
                  n_doubles_o,            &
                  wf%ao%n,                &
                  one,                    &
                  X,                      &
                  n_doubles_o,            &
                  wf%orbital_coefficients,&
                  wf%ao%n,                &
                  zero,                   &
                  wf%O_o,                 &
                  n_doubles_o)
!
      call mem%dealloc(X, n_doubles_o, wf%ao%n)
!
      call mem%alloc(X, n_doubles_v, wf%ao%n)
!
      call dgemm('T', 'N',                                      &
                  n_doubles_v,                                  &
                  wf%ao%n,                                      &
                  wf%ao%n,                                      &
                  one,                                          &
                  wf%orbital_coefficients_cc2(1, wf%n_o + 1),   &
                  wf%ao%n,                                      &
                  S,                                            &
                  wf%ao%n,                                      &
                  zero,                                         &
                  X,                                            &
                  n_doubles_v)
!
      call mem%dealloc(S, wf%ao%n, wf%ao%n)
!
      call dgemm('N', 'N',                               &
                  n_doubles_v,                           &
                  n_doubles_v,                           &
                  wf%ao%n,                               &
                  one,                                   &
                  X,                                     &
                  n_doubles_v,                           &
                  wf%orbital_coefficients(1, wf%n_o + 1),&
                  wf%ao%n,                               &
                  zero,                                  &
                  wf%O_v,                                &
                  n_doubles_v)
!
      call mem%dealloc(X, n_doubles_v, wf%ao%n)
!
   end subroutine construct_mlccsd_basis_transformation_matrix_mlccsd
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
!
      logical, intent(in), optional :: occupied_only
!
      logical  :: occupied_only_local
!
      real(dp), dimension(:,:), allocatable :: D
!
      integer, dimension(:), allocatable :: active_aos_ccsd, active_aos_cc2
!
      integer :: first_ao_ccsd, last_ao_ccsd, n_active_aos_ccsd
      integer :: first_ao_cc2, last_ao_cc2, n_active_aos_cc2
!
      integer :: i, mo_offset
!
      real(dp), parameter :: full_cd_threshold = 1.0d-4
!
      occupied_only_local = .false. 
!
      if (present(occupied_only)) occupied_only_local = occupied_only
!
!     Construct active occupied orbitals     
!
!     Determine CCSD ao list
!
!     Note that active atoms are the first atoms of the array of atoms in system. 
!     They are also ordered after the method they are treated with.
!
      call wf%ao%get_aos_in_subset('ccsd', first_ao_ccsd, last_ao_ccsd)
!
      n_active_aos_ccsd = last_ao_ccsd - first_ao_ccsd + 1
!
      call mem%alloc(active_aos_ccsd, n_active_aos_ccsd)
!
      do i = 1, n_active_aos_ccsd
!
         active_aos_ccsd(i) = first_ao_ccsd + i - 1
!
      enddo
!      
!     If we do both CC2 and CCS, then determine the CC2 ao list 
!     (from both CCSD and CC2 active atoms)
!
      if (wf%do_cc2 .and. wf%do_ccs) then
!
         call wf%ao%get_aos_in_subset('cc2', first_ao_cc2, last_ao_cc2)
!         
         n_active_aos_cc2 = last_ao_cc2 - first_ao_cc2 + 1
         call mem%alloc(active_aos_cc2, n_active_aos_ccsd + n_active_aos_cc2)
!
         do i = 1, n_active_aos_ccsd
!
            active_aos_cc2(i) = first_ao_ccsd + i - 1
!
         enddo
!
         do i = 1, n_active_aos_cc2 
!
            active_aos_cc2(i + n_active_aos_ccsd) = first_ao_cc2 + i - 1
!
         enddo
!
      endif
!
!     Set up active occupied density 
!
      call mem%alloc(D, wf%ao%n, wf%ao%n)
!
      call dgemm('N', 'T',                   &
                  wf%ao%n,                   &
                  wf%ao%n,                   &
                  wf%n_o,                    &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%ao%n,                   &
                  wf%orbital_coefficients,   &
                  wf%ao%n,                   &
                  zero,                      &
                  D,                         &
                  wf%ao%n)
!
!
!     Construct active occupied orbitals
!
      mo_offset = 0
!
      call wf%construct_orbital_block_by_density_cd(D, wf%n_ccsd_o, &
                                 wf%cholesky_orbital_threshold, mo_offset, active_aos_ccsd)
!
      mo_offset = mo_offset + wf%n_ccsd_o
!
      if (wf%do_cc2 .and. wf%do_ccs) then
!
         call wf%construct_orbital_block_by_density_cd(D, wf%n_cc2_o, &
                           wf%cholesky_orbital_threshold, mo_offset, active_aos_cc2)
!  
         mo_offset = mo_offset + wf%n_cc2_o
!
         call wf%construct_orbital_block_by_density_cd(D, wf%n_ccs_o, &
                           full_cd_threshold, mo_offset)
!
      elseif (wf%do_cc2 .and. .not. wf%do_ccs) then 
!
         call wf%construct_orbital_block_by_density_cd(D, wf%n_cc2_o, full_cd_threshold, mo_offset)
         wf%n_ccs_o = 0
!
      elseif (wf%do_ccs .and. .not. wf%do_cc2) then
!
         call wf%construct_orbital_block_by_density_cd(D, wf%n_ccs_o, full_cd_threshold, mo_offset)
         wf%n_cc2_o = 0
!
      endif  
!
!     Construct active virtual orbitals     
!
      if (.not. occupied_only_local) then
!
!        Set up virtual density     
!
         call dgemm('N', 'T',                                  &
                     wf%ao%n,                                  &
                     wf%ao%n,                                  &
                     wf%n_v,                                   &
                     one,                                      &
                     wf%orbital_coefficients(1, wf%n_o + 1),   &
                     wf%ao%n,                                  &
                     wf%orbital_coefficients(1, wf%n_o + 1),   &
                     wf%ao%n,                                  &
                     zero,                                     &
                     D,                                        &
                     wf%ao%n)
!
         mo_offset = wf%n_o
!
         call wf%construct_orbital_block_by_density_cd(D, wf%n_ccsd_v, &
                              wf%cholesky_orbital_threshold, wf%n_o, active_aos_ccsd)
!
         mo_offset = mo_offset + wf%n_ccsd_v
!
         if (wf%do_cc2 .and. wf%do_ccs) then
!
            call wf%construct_orbital_block_by_density_cd(D, wf%n_cc2_v, &
                              wf%cholesky_orbital_threshold, mo_offset, active_aos_cc2)
!
            call wf%ao%get_aos_in_subset('cc2', first_ao_cc2, last_ao_cc2)
!         
            n_active_aos_cc2 = last_ao_cc2 - first_ao_cc2 + 1
            call mem%dealloc(active_aos_cc2, n_active_aos_ccsd + n_active_aos_cc2)
!
            mo_offset = mo_offset + wf%n_cc2_v
!
            call wf%construct_orbital_block_by_density_cd(D, wf%n_ccs_v, &
                                                         full_cd_threshold, mo_offset)
!
         elseif (wf%do_cc2 .and. .not. wf%do_ccs) then 
!
            call wf%construct_orbital_block_by_density_cd(D, wf%n_cc2_v, &
                                                         full_cd_threshold, mo_offset)
            wf%n_ccs_v = 0
!
         elseif (wf%do_ccs .and. .not. wf%do_cc2) then
!
            call wf%construct_orbital_block_by_density_cd(D, wf%n_ccs_v, &
                                                         full_cd_threshold, mo_offset)
            wf%n_cc2_v = 0
!
         endif
!
      endif
!
      call mem%dealloc(active_aos_ccsd, n_active_aos_ccsd) 
      call mem%dealloc(D, wf%ao%n, wf%ao%n)        
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
      real(dp), dimension(:,:), allocatable :: D, S, PAO_coeff
!
      integer :: first_ao, last_ao, n_active_aos, rank
!
!     0. Determine active ao list
!
      call wf%ao%get_aos_in_subset('ccsd', first_ao, last_ao)
!
      n_active_aos = last_ao - first_ao + 1
!
!     1. Set up occupied density
!
      call mem%alloc(D, wf%ao%n, wf%ao%n)
!
      call dgemm('N', 'T',                   &
                  wf%ao%n,                   &
                  wf%ao%n,                   &
                  wf%n_o,                    &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%ao%n,                   &
                  wf%orbital_coefficients,   &
                  wf%ao%n,                   &
                  zero,                      &
                  D,                         &
                  wf%ao%n)
!
      if (wf%exists_frozen_fock_terms) then
!
         call daxpy(wf%ao%n**2, one, wf%frozen_CCT, 1, D, 1)
!
      endif
!
!     2. Construct PAOs for active atoms
!
      call mem%alloc(PAO_coeff, wf%ao%n, n_active_aos)
!
      call wf%project_atomic_orbitals(D, PAO_coeff, n_active_aos, first_ao)
!
!     3. Orthonormalize PAOs to get active virtual orbitals
!
      call mem%alloc(S, n_active_aos, n_active_aos)
!
      call wf%get_orbital_overlap(PAO_coeff, n_active_aos, S)
!
      call wf%lowdin_orthonormalization(PAO_coeff, S, n_active_aos, rank)
!
      call mem%dealloc(S, n_active_aos, n_active_aos)
!
      wf%n_ccsd_v = rank 
!
!     Set the active virtual orbital coefficients
!
      call dcopy(wf%n_ccsd_v*wf%ao%n, PAO_coeff, 1, wf%orbital_coefficients(1, wf%n_o + 1), 1)
!
      call mem%dealloc(PAO_coeff, wf%ao%n, n_active_aos)
!
      if (rank .lt. wf%n_v) then 
!
         if (wf%do_cc2 .and. wf%do_ccs) then
!
            call wf%ao%get_aos_in_subset('cc2', first_ao, last_ao)
!
            n_active_aos = last_ao - first_ao + 1
!
!           1. Set up occupied + virtual CCSD density 
!
            call dgemm('N', 'T',             &
                     wf%ao%n,                &
                     wf%ao%n,                &
                     wf%n_o + wf%n_ccsd_v,   &
                     one,                    &
                     wf%orbital_coefficients,&
                     wf%ao%n,                &
                     wf%orbital_coefficients,&
                     wf%ao%n,                &
                     zero,                   &
                     D,                      &
                     wf%ao%n)
!
            if (wf%exists_frozen_fock_terms) then
!  
               call daxpy(wf%ao%n**2, one, wf%frozen_CCT, 1, D, 1)
!  
            endif
!
!           2. Construct PAOs for active atoms
!
            call mem%alloc(PAO_coeff, wf%ao%n, n_active_aos)
!
            call wf%project_atomic_orbitals(D, PAO_coeff, n_active_aos, first_ao)
!
!           3. Orthonormalize PAOs to get active virtual orbitals
!
            call mem%alloc(S, n_active_aos, n_active_aos)
!
            call wf%get_orbital_overlap(PAO_coeff, n_active_aos, S)
!
            call wf%lowdin_orthonormalization(PAO_coeff, S, n_active_aos, rank)
!
            call mem%dealloc(S, n_active_aos, n_active_aos)
!
            wf%n_cc2_v = rank
!
!           Set the active virtual orbital coefficients
!
            call dcopy(wf%n_cc2_v*wf%ao%n, PAO_coeff, 1, &
                     wf%orbital_coefficients(1, wf%n_o + wf%n_ccsd_v + 1), 1)
!
            call mem%dealloc(PAO_coeff, wf%ao%n, n_active_aos)
!
            if (wf%n_ccsd_v + wf%n_cc2_v .lt. wf%n_v) then 
!
!              Set virtual inactive
!
!              4. Construct M = sum_p C_αp C_βp  for p = 1, n_o + n_cc2_v + n_ccsd_v
!
               call dgemm('N', 'T',                            &
                           wf%ao%n,                            &
                           wf%ao%n,                            &
                           wf%n_o + wf%n_ccsd_v + wf%n_cc2_v,  &
                           one,                                &
                           wf%orbital_coefficients,            &
                           wf%ao%n,                            &
                           wf%orbital_coefficients,            &
                           wf%ao%n,                            &
                           zero,                               &
                           D,                                  &
                           wf%ao%n)
!
               if (wf%exists_frozen_fock_terms) then
!
                  call daxpy(wf%ao%n**2, one, wf%frozen_CCT, 1, D, 1)
!
               endif
!
!              Construct PAOs for the remaining virtual orbitals
!
               call mem%alloc(PAO_coeff, wf%ao%n, wf%ao%n)
!
               call wf%project_atomic_orbitals(D, PAO_coeff, wf%ao%n)
!
!              5. Orthonormalize PAOs to get inactive virtual orbitals
!
               call mem%alloc(S, wf%ao%n, wf%ao%n)
!
               call wf%get_orbital_overlap(PAO_coeff, wf%ao%n, S)
!
               call wf%lowdin_orthonormalization(PAO_coeff, S, wf%ao%n, rank)
!
               call mem%dealloc(S, wf%ao%n, wf%ao%n)
!
               wf%n_ccs_v = rank 
!
!              6. Set inactive virtuals
!
               call dcopy(wf%n_ccs_v*wf%ao%n, PAO_coeff, 1, &
                           wf%orbital_coefficients(1, wf%n_o + wf%n_ccsd_v + wf%n_cc2_v + 1), 1)
!
               call mem%dealloc(PAO_coeff, wf%ao%n, wf%ao%n)
!
            endif
!
         else 
!
!           Set virtual inactive
!
!           4. Construct M = sum_p C_αp C_βp  for p = 1, n_o + n_ccsd_v
!
            call dgemm('N', 'T',             &
                     wf%ao%n,                &
                     wf%ao%n,                &
                     wf%n_o + wf%n_ccsd_v,   &
                     one,                    &
                     wf%orbital_coefficients,&
                     wf%ao%n,                &
                     wf%orbital_coefficients,&
                     wf%ao%n,                &
                     zero,                   &
                     D,                      &
                     wf%ao%n)
!
            if (wf%exists_frozen_fock_terms) then
!
               call daxpy(wf%ao%n**2, one, wf%frozen_CCT, 1, D, 1)
!
            endif
!
!           Construct PAOs for the remaining virtual orbitals
!
            call mem%alloc(PAO_coeff, wf%ao%n, wf%ao%n)
!
            call wf%project_atomic_orbitals(D, PAO_coeff, wf%ao%n)
!
!           5. Orthonormalize PAOs to get inactive virtual orbitals
!
            call mem%alloc(S, wf%ao%n, wf%ao%n)
!
            call wf%get_orbital_overlap(PAO_coeff, wf%ao%n, S)
!
            call wf%lowdin_orthonormalization(PAO_coeff, S, wf%ao%n, rank)
!
            call mem%dealloc(S, wf%ao%n, wf%ao%n)
!
            if (wf%do_cc2) then
!
               wf%n_cc2_v = rank 
!
            elseif (wf%do_ccs) then
!
               wf%n_ccs_v = rank 
!
            endif
!
!           6. Set inactive virtuals
!
            call dcopy((rank)*wf%ao%n, PAO_coeff, 1, &
                           wf%orbital_coefficients(1, wf%n_o + wf%n_ccsd_v + 1), 1)
!
            call mem%dealloc(PAO_coeff, wf%ao%n, wf%ao%n)
!
         endif
!
      endif
!
      call mem%dealloc(D, wf%ao%n, wf%ao%n)
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
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(out) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(out) :: T_v
!
      integer :: k
!
      logical :: set_to_zero
!
      integer :: n_cnto_states
!
      character(len=200) :: r_or_l
!
      real(dp), dimension(:,:,:), allocatable   :: R_ai
      real(dp), dimension(:,:), allocatable     :: R_aibj
      real(dp), dimension(:,:,:,:), allocatable :: R_aibj_sq
!
      r_or_l = 'right'
!
      if (input%is_keyword_present('left eigenvectors', 'solver cc es')) r_or_l = 'left'
!
      if(wf%do_cc2 .and. wf%do_ccs) then
!
         call output%error_msg('CNTOs can not be constructed for three-level MLCCSD.' //&
                               ' Try using cnto-approx in stead.')
!
      endif
!
!     Run CC2 calculation
!
      n_cnto_states = size(wf%cnto_states)
!
      call mem%alloc(R_ai, wf%n_v, wf%n_o, n_cnto_states)
      call mem%alloc(R_aibj, wf%n_t1*(wf%n_t1+1)/2, n_cnto_states)
!
      call wf%cc2_calculation_for_cntos(r_or_l, n_cnto_states, &
                                          R_ai, R_aibj, wf%cnto_states)
!
      call mem%alloc(R_aibj_sq, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      set_to_zero = .true.
!
      do k = 1, n_cnto_states
!
         call squareup(R_aibj(:,k), R_aibj_sq, wf%n_t1)
!
!        Add contribution to M and N
!
         call wf%construct_M_and_N_cnto(R_ai(:,:,k), R_aibj_sq, T_o, T_v, set_to_zero)
!
         set_to_zero = .false.
!
      enddo
!
      call mem%dealloc(R_ai, wf%n_v, wf%n_o, n_cnto_states)
      call mem%dealloc(R_aibj, wf%n_t1*(wf%n_t1+1)/2, n_cnto_states)
      call mem%dealloc(R_aibj_sq, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%diagonalize_M_and_N(T_o, T_v)
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
!
      use diis_cc_gs_class, only: diis_cc_gs
      use davidson_cc_es_class, only: davidson_cc_es
      use cc2_class, only: cc2
!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      character(len=200), intent(in) :: transformation
!
      integer, intent(in) :: n_cnto_states
!
      real(dp), dimension(wf%n_v, wf%n_o, n_cnto_states), intent(out) :: R_ai
      real(dp), dimension(wf%n_t1*(wf%n_t1+1)/2, n_cnto_states), intent(out) :: R_aibj
!
      integer, dimension(n_cnto_states), intent(in) :: cnto_states
!
      type(cc2), allocatable :: cc2_wf
!
      type(diis_cc_gs), allocatable :: cc_gs_solver_diis
      type(davidson_cc_es), allocatable :: cc_es_solver_davidson
!
      type(timings) :: timer_gs, timer_es
!
      real(dp), dimension(:), allocatable :: R
!  
      integer :: n, n_es
!
      call output%printf('m', 'Running CC2 calculation for CNTOs.', fs='(/t3,a)')
!
      if (.not. input%is_keyword_present('print cc2 calculation','mlcc')) &
         call output%mute()
!
!     Run CC2 calculation
!
      allocate(cc2::cc2_wf)
      call cc2_wf%initialize(wf)
!
      call cc2_wf%mo_preparations()
!
      cc2_wf%eri = t1_eri_tool(wf%eri)
      call cc2_wf%eri%initialize()
      call cc2_wf%eri%copy_from_t1(wf%eri)
!
!     1. Ground state
!
      timer_gs = timings('Ground state CC2 calculation for CNTOs')
      call timer_gs%turn_on()
!
      cc_gs_solver_diis = diis_cc_gs(cc2_wf, restart=.false.)
      call cc_gs_solver_diis%run(cc2_wf)
      call cc_gs_solver_diis%cleanup(cc2_wf)
      call timer_gs%turn_off()
!
!     Excited states
!
      timer_es = timings('Excited state CC2 calculation for CNTOs')
      call timer_es%turn_on()
!
      call cc2_wf%construct_fock('es')
!
      cc_es_solver_davidson = davidson_cc_es(transformation, cc2_wf, restart=.false.)
      call cc_es_solver_davidson%run(cc2_wf)
      call cc_es_solver_davidson%cleanup(cc2_wf)
!
      call timer_es%turn_off()
!
!     Transfer information to mlcc wavefunction
!
!     2. Excitation vectors
!
      n_es = cc2_wf%n_singlet_states
!
      if(n_es .lt. n_cnto_states) call output%error_msg('Requested too many CNTO states')
!
!     Return only the excitation vectors in cnto_states
!
      call mem%alloc(R, cc2_wf%n_t1 + cc2_wf%n_t2)
!
      do n = 1, n_cnto_states
!
         if(n_es .lt. cnto_states(n)) call output%error_msg('Requested non-existent CNTO state')
!
         call cc2_wf%read_excited_state(R,               &
                                        cnto_states(n),  &
                                        cnto_states(n),  &
                                        transformation)   
!
         call dcopy(cc2_wf%n_t1, R, 1, R_ai(1, 1,n), 1)
         call dcopy(cc2_wf%n_t2, R(cc2_wf%n_t1 + 1), 1, R_aibj(1,n), 1)
!
      enddo  
!
      call mem%dealloc(R, cc2_wf%n_t1 + cc2_wf%n_t2)
!
!     Cleanup and print
!
      call cc2_wf%cleanup() 
!
      if (.not. input%is_keyword_present('print cc2 calculation','mlcc')) &
         call output%unmute()
!
      call output%printf('m', '- Summary of CC2 calculation for CNTOs:',fs='(/t3,a)')
!
      call output%printf('m', 'Wall time for CC2 ground calculation (sec):   (f20.2)', &
                         reals=[timer_gs%get_elapsed_time('wall')], fs='(/t6,a)')
      call output%printf('m', 'CPU time for CC2 ground calculation (sec):    (f20.2)', &
                         reals=[timer_gs%get_elapsed_time('cpu')], fs='(t6,a)')
!
      call output%printf('m', 'Wall time for CC2 excited calculation (sec):  (f20.2)', &
                         reals=[timer_es%get_elapsed_time('wall')], fs='(/t6,a)')
      call output%printf('m', 'CPU time for CC2 excited calculation (sec):   (f20.2)', &
                         reals=[timer_es%get_elapsed_time('cpu')], fs='(t6,a)')
!
   end subroutine cc2_calculation_for_cntos_mlccsd
!
!
   module subroutine construct_semicanonical_mlcc_orbitals_mlccsd(wf)
!!
!!    Construct semicanonical orbitals
!!    Written by Sarai D. Folkestad
!!
!!    Wrapper to construct orbitals that block diagonalizes the fock matrix
!!    for the different orbitals
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      integer :: n_levels
!
      n_levels = 3
!
      call wf%construct_block_diagonal_fock_orbitals(n_levels,                    &
                                          [wf%n_ccsd_o, wf%n_cc2_o, wf%n_ccs_o],  &
                                          [wf%n_ccsd_v, wf%n_cc2_v, wf%n_ccs_v],  &
                                          wf%orbital_coefficients,                &
                                          wf%orbital_energies) 
!
   end subroutine construct_semicanonical_mlcc_orbitals_mlccsd
!
!
   module subroutine construct_orbitals_cc2_mlccsd(wf)
!!
!!    Construct orbitals cc2
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      integer :: n_levels
!
      n_levels = 2
!
      call wf%initialize_t1()
      call zero_array(wf%t1, wf%n_t1)
      call wf%eri%set_t1_to_mo()
!
      call wf%construct_fock()
      call wf%destruct_t1()

      call dcopy(wf%n_mo*wf%ao%n, wf%orbital_coefficients, 1, wf%orbital_coefficients_cc2, 1)
!
      call wf%construct_block_diagonal_fock_orbitals(n_levels,                   &
                                       [wf%n_ccsd_o + wf%n_cc2_o, wf%n_ccs_o],   &
                                       [wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccs_v],   &
                                       wf%orbital_coefficients_cc2,              &
                                       wf%orbital_energies_cc2) 
!
      call wf%initialize_O_o()
      call wf%initialize_O_v()
!
      call wf%construct_mlccsd_basis_transformation_matrix()
!
   end subroutine construct_orbitals_cc2_mlccsd
!
!
end submodule orbitals_mlccsd
