submodule (mlccsd_class) orbital_partitioning
!
!!
!!    Orbital partitioning submodule (MLCCSD) 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2017
!!
!!    Contains the following family of procedures of the MLCCSD class:
!!
!!
!
   implicit none 
!
   logical :: debug   = .false.
   logical :: timings = .true.
!
!
contains
!
!
   module subroutine cholesky_localization_drv_mlccsd(wf)
!!
!!    Cholesky orbital localization. driver,
!!    Written by Sarai D. Folkestad, June 2017
!!

      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(:,:), allocatable :: orbitals
      real(dp), dimension(:,:), allocatable :: orbital_energies
      real(dp), dimension(:,:), allocatable :: S
      real(dp), dimension(:,:), allocatable :: S_packed
      real(dp), dimension(:,:), allocatable :: Y_1
      real(dp), dimension(:,:), allocatable :: Y_2
!
      integer(i15) :: unit_overlap
!  
      integer(i15) :: n_nuclei
      integer(i15) :: active_space = 0, n_CC2_atoms = 0, n_CCSD_atoms = 0
      integer(i15) :: n_active_orbitals_o = 0, n_active_orbitals_v = 0 
      integer(i15) :: offset_o = 0, offset_v = 0
      integer(i15) :: n_active_aos = 0
!
      integer(i15), dimension(:,:), allocatable :: ao_center_info, n_ao_on_center
      integer(i15), dimension(:,:), allocatable :: active_atoms, active_atoms_CC2, active_atoms_CCSD
      integer(i15), dimension(:,:), allocatable :: active_ao_index_list
      integer(i15) :: n_vectors_o, n_vectors_v
!
!     Timing variables
!
      real(dp) :: start_chol_deco = 0, end_chol_deco = 0
!
!     IO-variables
!
      logical      :: file_exists     
      integer(i15) :: unit_cholesky_decomp = 0, ioerror = 0
!
!     Indices
!
      integer(i15) :: i = 0, j = 0, ij = 0
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: density_o, density_v
      real(dp), dimension(:,:), allocatable :: ao_fock
      real(dp), dimension(:,:), allocatable :: C_cc2
!
!     Start timings
!
      call cpu_time(start_chol_deco)
!
!     :::::::::::::::::::::::::::::::::::::
!     -::-  Prepare for decomposition  -::-
!     :::::::::::::::::::::::::::::::::::::
!
!     Prints
!
      write(unit_output,'(/t3,a/)')    ':: Cholesky decomposition '
      flush(unit_output)
!
!     Get center info
!
      call read_atom_info(n_nuclei, wf%n_ao)
!
      call allocator_int(n_ao_on_center, n_nuclei, 1)      
      call allocator_int(ao_center_info, wf%n_ao, 2)
!
      call read_center_info(n_nuclei, wf%n_ao, n_ao_on_center, ao_center_info)
!
!     Check that cholesky.inp exists
!
      inquire(file='cholesky.inp', exist=file_exists)
      if (.not. file_exists) then
         write(unit_output,*) 'WARNING: Input file for cholesky decomposition is not found.'
         stop
      endif
!
!     Open cholesky.inp
!
      call generate_unit_identifier(unit_cholesky_decomp)
      open(unit=unit_cholesky_decomp, file='cholesky.inp', status='unknown', form='formatted', iostat=ioerror)
      if (ioerror .ne. 0) write(unit_output,*)'WARNING: Error while opening cholesky.inp'
      rewind(unit_cholesky_decomp)
!
!     Start timings
!
      call cpu_time(start_chol_deco)
!
!     Construct AO fock matrix
!
      call allocator(density_o, wf%n_ao, wf%n_ao)
      call allocator(density_v, wf%n_ao, wf%n_ao)
      density_o = zero
      density_v = zero
!
      call wf%construct_density_matrices(density_o, density_v, wf%mo_coef, wf%n_o, wf%n_v)
!
      call allocator(ao_fock, wf%n_ao, wf%n_ao)
      call wf%construct_ao_fock(ao_fock)
!
      call deallocator(density_o, wf%n_ao, wf%n_ao)
      call deallocator(density_v, wf%n_ao, wf%n_ao)
!
!     Test for CCS region
!    
      if (wf%mlcc_settings%CCS) then
!
!        Test for CC2 region
!
         if (wf%mlcc_settings%CC2) then
            write(unit_output,*)'CCSD/CC2/CCS wavefunction requested'
            flush(unit_output)
            call wf%cholesky_localization_CCSD_CC2_CCS(ao_center_info, n_ao_on_center,&
                                                       ao_fock, n_nuclei, unit_cholesky_decomp)
!
         else ! No CC2 region
!
            write(unit_output,*)'CCS/CCSD wavefunction requested' 
            flush(unit_output)
            write(unit_output,*)'Partitioning not tested'
            stop
!
            call wf%cholesky_localization_CCSD_CCS(ao_center_info, n_ao_on_center,&
                                                       ao_fock, n_nuclei, unit_cholesky_decomp)
!
         endif
!
      else ! No CCS region
!
         write(unit_output,*)'CC2/CCSD wavefunction requested' 
         flush(unit_output)
!
         call wf%cholesky_localization_CCSD_CC2(ao_center_info, n_ao_on_center,&
                                                   ao_fock, n_nuclei, unit_cholesky_decomp)
!
      endif
!
!     :: Check orthogonality of new orbitals ::
!
!     Read overlap matrix
!
      call generate_unit_identifier(unit_overlap)
      open(unit=unit_overlap, file='MLCC_OVERLAP', status='unknown', form='formatted', iostat=ioerror)
      if (ioerror .ne. 0) write(unit_output,*)'WARNING: Error while opening MLCC_OVERLAP'
      rewind(unit_overlap)
!
      call allocator(S_packed, wf%n_ao*(wf%n_ao+1)/2, 1)
      read(unit_overlap,*)S_packed
      close(unit_overlap)
!
      call allocator(S, wf%n_ao, wf%n_ao)
!
      call squareup(S_packed, S, wf%n_ao)
!
      call deallocator(S_packed, wf%n_ao*(wf%n_ao+1)/2, 1)
!
      call allocator(Y_1, wf%n_mo, wf%n_ao)
      call dgemm('T', 'N',    &
                  wf%n_mo,    &
                  wf%n_ao,    &
                  wf%n_ao,    &
                  one, &
                  wf%mo_coef, &
                  wf%n_ao,&
                  S, &
                  wf%n_ao, &
                  zero, &
                  Y_1, &
                  wf%n_mo)
!
      call deallocator(S, wf%n_ao, wf%n_ao)
      call allocator(Y_2, wf%n_mo, wf%n_mo)
!
      call dgemm('N', 'N',    &
                  wf%n_mo,    &
                  wf%n_mo,    &
                  wf%n_ao,    &
                  one, &
                  Y_1, &
                  wf%n_mo,&
                  wf%mo_coef, &
                  wf%n_ao, &
                  zero, &
                  Y_2, &
                  wf%n_mo)
!
      call deallocator(Y_1, wf%n_mo, wf%n_ao)
!
      if (.not. (check_orthogonality(Y_2, wf%n_mo, wf%n_mo))) then
         write(unit_output,*)'New orbitals not orthogonal'
         stop
      endif
!
      call deallocator(Y_2, wf%n_mo, wf%n_mo)
!
!
!     Deallocate ao-fock
!
      call deallocator(ao_fock, wf%n_ao, wf%n_ao)
      call deallocator_int(n_ao_on_center, n_nuclei, 1)      
      call deallocator_int(ao_center_info, wf%n_ao, 2)
!
!     Print decomposition info
!
      write(unit_output, '(/t3,a15/)')'Active space:  '
      write(unit_output,'(t3,a45, i3)') 'Number of CCSD active occupied orbitals:  ', wf%n_CCSD_o
      write(unit_output,'(t3,a45, i3/)')'Number of CCSD active virtual orbitals:   ', wf%n_CCSD_v
      write(unit_output,'(t3,a45, i3)') 'Number of CC2 active occupied orbitals:   ', wf%n_CC2_o
      write(unit_output,'(t3,a45, i3/)')'Number of CC2 active virtual orbitals:    ', wf%n_CC2_v                           
!
      write(unit_output, '(/t3,a15/)')'Inactive space:  '
      write(unit_output,'(t3,a45, i3)') 'Number of CCS inactive occupied orbitals: ', wf%n_CCS_o
      write(unit_output,'(t3,a45, i3/)')'Number of CCS inactive virtual orbitals:  ', wf%n_CCS_v
!
!     Print timings
!
      call cpu_time(end_chol_deco)
      if (timings) write(unit_output,'(/t3,a27,f14.8/)') 'Total time (seconds):', end_chol_deco - start_chol_deco
      flush(unit_output)
!
   end subroutine cholesky_localization_drv_mlccsd
!
   module subroutine cholesky_localization_CCSD_CC2_CCS_mlccsd(wf, ao_center_info, n_ao_on_center,&
                                                       ao_fock, n_nuclei, unit_cholesky_decomp)
!!
!!    Cholesky orbital localization. driver,
!!    Written by Sarai D. Folkestad, June 2017
!!
!!    Driver for Cholesky density decomposition.  
!!
      implicit none
!
!     Input arguments
!
      class(mlccsd) :: wf
!
      integer(i15), dimension(wf%n_ao, 2)  :: ao_center_info
      integer(i15), dimension(n_nuclei, 1) :: n_ao_on_center
!
      real(dp), dimension(:,:) :: ao_fock
!
      integer(i15) :: n_nuclei
      integer(i15) :: unit_cholesky_decomp
      integer(i15) :: i, j, ij
!
!     Local routine variables
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: orbitals
      real(dp), dimension(:,:), allocatable :: orbital_energies
      real(dp), dimension(:,:), allocatable :: density_o
      real(dp), dimension(:,:), allocatable :: density_v
      real(dp), dimension(:,:), allocatable :: C_cc2
!
      integer(i15), dimension(:,:), allocatable :: active_atoms, active_atoms_CC2, active_atoms_CCSD
      integer(i15), dimension(:,:), allocatable :: active_ao_index_list
!
!     Active space variables
!
      integer(i15) :: n_active_aos
      integer(i15) :: n_active_orbitals_o
      integer(i15) :: n_active_orbitals_v
      integer(i15) :: n_CC2_atoms, n_CCSD_atoms
      integer(i15) :: n_vectors_o, n_vectors_v
      integer(i15) :: offset_o, offset_v, offset
!
!     :::::::::::::::::::::::::::::::::::::
!     -::-  Prepare for decomposition  -::-
!     :::::::::::::::::::::::::::::::::::::
!
      call allocator(orbitals, wf%n_ao, wf%n_mo)
      call allocator(orbital_energies, wf%n_mo, 1)
!
!     ::::::::::::::::::::::::::::::::
!     -::-    CC2/CCS partition   -::-
!     ::::::::::::::::::::::::::::::::
!
!     :: Construct canonical occupied and vacant density matrices ::    
!
      call allocator(density_o, wf%n_ao, wf%n_ao)
      call allocator(density_v, wf%n_ao, wf%n_ao)
      density_o = zero
      density_v = zero
!
      call wf%construct_density_matrices(density_o, density_v, wf%mo_coef, wf%n_o, wf%n_v)
!
      offset_o = 1
      offset_v = 1 + wf%n_o
!  
!     Get CC2/CCSD-active atoms
!  
      n_CC2_atoms  =  get_number_of_active_atoms(unit_cholesky_decomp, 'CC2  ')
      n_CCSD_atoms =  get_number_of_active_atoms(unit_cholesky_decomp, 'CCSD ')
!  
!     Allocate active atoms list
!  
      call allocator_int(active_atoms_CC2, n_CC2_atoms, 1)
      call allocator_int(active_atoms_CCSD, n_CCSD_atoms, 1)
      call allocator_int(active_atoms, n_CCSD_atoms + n_CC2_atoms, 1)
!  
!     Get list of active atoms
! 
      call get_active_atoms(unit_cholesky_decomp, active_atoms_CC2,  n_CC2_atoms,  'CC2  ')
      call get_active_atoms(unit_cholesky_decomp, active_atoms_CCSD, n_CCSD_atoms, 'CCSD ')
!
      do i = 1, n_CC2_atoms
         active_atoms(i, 1) = active_atoms_CC2(i, 1)
      enddo
      do i = 1, n_CCSD_atoms
         active_atoms(n_CC2_atoms + i, 1) = active_atoms_CCSD(i,1)
      enddo
!
      call deallocator_int(active_atoms_CC2, n_CC2_atoms, 1)
      call deallocator_int(active_atoms_CCSD, n_CCSD_atoms, 1)
! 
!     Sanity check on active atoms
! 
       if ((n_CC2_atoms + n_CCSD_atoms) .gt. n_nuclei) then
         write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
         stop
      endif
! 
      do i = 1, n_CC2_atoms + n_CCSD_atoms
         if (active_atoms(i,1) .gt. n_nuclei) then
            write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.',i, active_atoms(i,1)
            stop
         endif
      enddo
!  
!     :: CC2 localized Cholesky orbitals ::
!  
      n_active_aos = 0
!
      do i = 1, n_CC2_atoms + n_CCSD_atoms
         n_active_aos = n_active_aos + n_ao_on_center(active_atoms(i,1),1)
      enddo
!  
!     Construct active_ao_index_list
!  
      call allocator_int(active_ao_index_list, n_active_aos, 1)
      call construct_active_ao_index_list(active_ao_index_list, n_active_aos, active_atoms, &
                                               n_CC2_atoms + n_CCSD_atoms, ao_center_info, wf%n_ao)
!  
!     Occupied part
!
      n_vectors_o = 0  
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                           density_o, n_vectors_o,&
                                          .true., n_active_aos, active_ao_index_list)
!  
!     Virtual part
!  
      n_vectors_v = 0 
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                           density_v, n_vectors_v,&
                                           .true., n_active_aos, active_ao_index_list)
!
!
!     Save CC2 space information
!
      wf%n_CC2_o = n_vectors_o
      wf%n_CC2_v = n_vectors_v
!  
!     Calculate new offset         
!  
      offset_o = offset_o + n_vectors_o
      offset_v = offset_v + n_vectors_v
!  
      call deallocator_int(active_atoms, n_CC2_atoms + n_CCSD_atoms, 1)
      call deallocator_int(active_ao_index_list, n_active_aos, 1)
!     
!
!     :: CCS localized Cholesky orbitals  ::
!
!     Occupied part
!
      n_vectors_o = 0 
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                        density_o, n_vectors_o,&
                                        .false., n_active_aos)
!
!     Virtual part
!
      n_vectors_v = 0 
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                        density_v, n_vectors_v,&
                                        .false., n_active_aos)
!
!     Save inactive space information
!
      wf%n_CCS_o = n_vectors_o
      wf%n_CCS_v = n_vectors_v
!
      wf%first_CCS_o = offset_o
      wf%first_CCS_v = offset_v - wf%n_o
!
!     Totsal number of active orbitals
!
      n_active_orbitals_o = wf%n_o - wf%n_CCS_o
      n_active_orbitals_v = wf%n_v - wf%n_CCS_v
!
      wf%mo_coef_cc2_ccs       = orbitals
      wf%fock_diagonal_cc2_ccs = orbital_energies
!
!     Deallocations
!
      call deallocator(density_v, wf%n_ao, wf%n_ao)   
      call deallocator(density_o, wf%n_ao, wf%n_ao)
!
!     ::::::::::::::::::::::::::::::::
!     -::-   CCSD/CC2 partition   -::-
!     ::::::::::::::::::::::::::::::::
!
      offset_o = 1
      offset_v = 1 + wf%n_o
!
!     Construct CC2 density matrix, occupied and virtual
!
      call allocator(density_v, wf%n_ao, wf%n_ao)   
      call allocator(density_o, wf%n_ao, wf%n_ao)
!
!     Construct CC2 density
!
      call allocator(C_cc2, wf%n_ao*(wf%n_CC2_o + wf%n_CC2_v), 1) 
!
      do i = 1, wf%n_ao
         do j = 1, wf%n_CC2_o
            ij =  index_two(i, j, wf%n_ao)
            C_cc2(ij, 1) = wf%mo_coef_cc2_ccs(i, j)
         enddo
      enddo
!
      do i = 1, wf%n_ao
         do j = 1, wf%n_CC2_v
            ij =  index_two(i, wf%n_CC2_o + j, wf%n_ao)
            C_cc2(ij, 1) = wf%mo_coef_cc2_ccs(i, wf%n_o + j)
         enddo
      enddo
!
      call wf%construct_density_matrices(density_o, density_v, C_cc2, &
                                    wf%n_CC2_o, wf%n_CC2_v)
!  
      call deallocator(C_cc2, wf%n_ao*(wf%n_CC2_o + wf%n_CC2_v), 1) 
!
      n_CCSD_atoms =  get_number_of_active_atoms(unit_cholesky_decomp, 'CCSD ')
!  
!     Allocate active atoms list
!
      call allocator_int(active_atoms, n_CCSD_atoms, 1)
!  
!     Get list of active atoms
! 
      call get_active_atoms(unit_cholesky_decomp, active_atoms,  n_CCSD_atoms, 'CCSD ')
!  
!     Construct active_ao_index_list
!  
      n_active_aos = 0
!
      do i = 1, n_CCSD_atoms
         n_active_aos = n_active_aos + n_ao_on_center(active_atoms(i,1),1)
      enddo
!
      call allocator_int(active_ao_index_list, n_active_aos, 1)
      active_ao_index_list = 0
!
      call construct_active_ao_index_list(active_ao_index_list, n_active_aos, active_atoms, &
                                               n_CCSD_atoms, ao_center_info, wf%n_ao)
!
      call deallocator_int(active_atoms, n_CCSD_atoms, 1)
!  
!     Occupied part
!  
      n_vectors_o = 0 
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                           density_o, n_vectors_o,&
                                          .true., n_active_aos, active_ao_index_list)
!  
!     Virtual part
!
      n_vectors_v = 0 
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                           density_v, n_vectors_v,&
                                           .true., n_active_aos, active_ao_index_list)
!
      call deallocator_int(active_ao_index_list, n_active_aos, 1)
!
!  
!     Save active space information
!  
      wf%n_CCSD_o = n_vectors_o
      wf%n_CCSD_v = n_vectors_v
!
      wf%first_CCSD_o = offset_o
      wf%first_CCSD_v = offset_v - wf%n_o
!  
!     Calculate new offset         
!  
      offset_o = offset_o + n_vectors_o
      offset_v = offset_v + n_vectors_v  
!
!     :: CC2  localized Cholesky orbitals  ::
!  
!     Occupied part
!  
      n_vectors_o = 0 
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                           density_o, n_vectors_o,&
                                          .false., n_active_aos)
!  
!     Virtual part
!  
      n_vectors_v = 0 
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                           density_v, n_vectors_v,&
                                           .false., n_active_aos)
!
!     Save inactive space information
!
      wf%n_CC2_o = n_vectors_o
      wf%n_CC2_v = n_vectors_v
!
      wf%first_CC2_o = offset_o
      wf%first_CC2_v = offset_v - wf%n_o
!  
!     Calculate new offset         
!  
      offset_o = offset_o + n_vectors_o
      offset_v = offset_v + n_vectors_v
!
      call deallocator(density_v, wf%n_ao, wf%n_ao)   
      call deallocator(density_o, wf%n_ao, wf%n_ao)
!
      wf%mo_coef       = orbitals
      wf%fock_diagonal = orbital_energies
!
      call deallocator (orbitals, wf%n_ao, wf%n_mo)
      call deallocator (orbital_energies, wf%n_mo, 1)
!
   end subroutine cholesky_localization_CCSD_CC2_CCS_mlccsd
!
!
   module subroutine cholesky_localization_CCSD_CCS_mlccsd(wf, ao_center_info, n_ao_on_center,&
                                                       ao_fock, n_nuclei, unit_cholesky_decomp)
!!
!!    Cholesky orbital localization. driver,
!!    Written by Sarai D. Folkestad, June 2017
!!
!!    Driver for Cholesky density decomposition.  
!!
!!
      implicit none
!
!     Input arguments
!
      class(mlccsd) :: wf
!
      integer(i15), dimension(wf%n_ao, 2)  :: ao_center_info
      integer(i15), dimension(n_nuclei, 1) :: n_ao_on_center
!
      real(dp), dimension(:,:) :: ao_fock
!
      integer(i15) :: n_nuclei
      integer(i15) :: unit_cholesky_decomp
!
!     Local routine variables
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: orbitals
      real(dp), dimension(:,:), allocatable :: orbital_energies
      real(dp), dimension(:,:), allocatable :: density_o
      real(dp), dimension(:,:), allocatable :: density_v
!
      integer(i15), dimension(:,:), allocatable :: active_atoms
      integer(i15), dimension(:,:), allocatable :: active_ao_index_list
!
!     Looping variables
!
      integer(i15) :: i = 0, j = 0, ij = 0
!
!
!
      integer(i15) :: n_active_aos
      integer(i15) :: n_active_orbitals_o
      integer(i15) :: n_active_orbitals_v
      integer(i15) :: n_CCSD_atoms
      integer(i15) :: n_vectors_o, n_vectors_v
      integer(i15) :: offset_o, offset_v
!
!     :::::::::::::::::::::::::::::
!     -::-  CCS/CCSD partition -::-
!     :::::::::::::::::::::::::::::
!
!
      call allocator(orbitals, wf%n_ao, wf%n_mo)
      call allocator(orbital_energies, wf%n_mo, 1)
!
!     :: Construct canonical occupied and vacant density matrices ::    
!
      call allocator(density_o, wf%n_ao, wf%n_ao)
      call allocator(density_v, wf%n_ao, wf%n_ao)
      density_o = zero
      density_v = zero
!
      call wf%construct_density_matrices(density_o, density_v, wf%mo_coef, wf%n_o, wf%n_v)
!
!     No CC2 region
!
      wf%n_CC2_o = 0
      wf%n_CC2_v = 0
!
      offset_o = 1
      offset_v = 1 + wf%n_o
!
!  
!     Get CC2/CCSD-active atoms
!
      n_CCSD_atoms =  get_number_of_active_atoms(unit_cholesky_decomp, 'CCSD ')  
!
!     Allocate active atoms list
!  
      call allocator_int(active_atoms, n_CCSD_atoms, 1)
!  
!     Get list of active atoms
!  
      call get_active_atoms(unit_cholesky_decomp, active_atoms,  n_CCSD_atoms, 'CCSD  ')
!  
!     Sanity check on active atoms
!  
      if (n_CCSD_atoms .gt. n_nuclei) then
         write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
         stop
      endif
!  
      do i = 1, n_CCSD_atoms
         if (active_atoms(i,1) .gt. n_nuclei) then
            write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
            stop
         endif
      enddo
!  
!     :: Constructing active (CCSD) localized Cholesky orbitals ::
!  
      n_active_aos = 0
      do i = 1, n_CCSD_atoms
         n_active_aos = n_active_aos + n_ao_on_center(active_atoms(i,1),1)
      enddo
!  
!     Construct active_ao_index_list
!  
      call allocator_int(active_ao_index_list, n_active_aos, 1)
      call construct_active_ao_index_list(active_ao_index_list, n_active_aos, active_atoms, &
                                               n_CCSD_atoms, ao_center_info, wf%n_ao)
!  
!     Occupied part
!  
      n_vectors_o = 0 
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                           density_o, n_vectors_o,&
                                          .true., n_active_aos, active_ao_index_list)
!  
!     Virtual part
!  
      n_vectors_v = 0 
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                           density_v, n_vectors_v,&
                                           .true., n_active_aos, active_ao_index_list)
!  
!     Save active space information
!  
      wf%n_CCSD_o = n_vectors_o
      wf%n_CCSD_v = n_vectors_v
!
      wf%first_CCSD_o = offset_o
      wf%first_CCSD_v = offset_v - wf%n_o
!  
!     Calculate new offset         
!  
      offset_o = offset_o + n_vectors_o
      offset_v = offset_v + n_vectors_v
!  
      call deallocator_int(active_atoms, n_CCSD_atoms, 1)
      call deallocator_int(active_ao_index_list, n_active_aos, 1)
!
!     :: CCS  localized Cholesky orbitals  ::
!
!     Occupied part
!
      n_vectors_o = 0 
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                        density_o, n_vectors_o,&
                                        .false., n_active_aos)
!
!     Virtual part
!
      n_vectors_v = 0 
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                        density_v, n_vectors_v,&
                                        .false., n_active_aos)
!
!     Save inactive space information
!
      wf%n_CCS_o = n_vectors_o
      wf%n_CCS_v = n_vectors_v
!
      wf%first_CCS_o = offset_o
      wf%first_CCS_v = offset_v - wf%n_o
!
      wf%mo_coef       = orbitals
      wf%fock_diagonal = orbital_energies
!
!     Deallocations
!
      call deallocator(density_v, wf%n_ao, wf%n_ao)   
      call deallocator(density_o, wf%n_ao, wf%n_ao)
!
      call deallocator (orbitals, wf%n_ao, wf%n_mo)
      call deallocator (orbital_energies, wf%n_mo, 1) 
!
   end subroutine cholesky_localization_CCSD_CCS_mlccsd
!
!
   module subroutine cholesky_localization_CCSD_CC2_mlccsd(wf, ao_center_info, n_ao_on_center,&
                                                       ao_fock, n_nuclei, unit_cholesky_decomp)
!!
!!    Cholesky orbital localization. driver,
!!    Written by Sarai D. Folkestad, June 2017
!!
!!    Driver for Cholesky density decomposition.  
!!
!!
      implicit none
!
!     Input arguments
!
      class(mlccsd) :: wf
!
      integer(i15), dimension(wf%n_ao, 2)  :: ao_center_info
      integer(i15), dimension(n_nuclei, 1) :: n_ao_on_center
!
      real(dp), dimension(:,:) :: ao_fock
!
      integer(i15) :: n_nuclei
      integer(i15) :: unit_cholesky_decomp
!
!     Local routine variables
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: orbitals
      real(dp), dimension(:,:), allocatable :: orbital_energies
      real(dp), dimension(:,:), allocatable :: density_o
      real(dp), dimension(:,:), allocatable :: density_v
!
      integer(i15), dimension(:,:), allocatable :: active_atoms
      integer(i15), dimension(:,:), allocatable :: active_ao_index_list
!
!     Looping variables
!
      integer(i15) :: i = 0, j = 0, ij = 0
!
!
!
      integer(i15) :: n_active_aos
      integer(i15) :: n_active_orbitals_o
      integer(i15) :: n_active_orbitals_v
      integer(i15) :: n_CCSD_atoms
      integer(i15) :: n_vectors_o, n_vectors_v
      integer(i15) :: offset_o, offset_v
!
      call allocator(orbitals, wf%n_ao, wf%n_mo)
      call allocator(orbital_energies, wf%n_mo, 1)
!      
      do i = 1, wf%n_ao
         do j = 1, wf%n_mo
            ij = index_two(i, j, wf%n_ao)
            wf%mo_coef_cc2_ccs (i,j)= wf%mo_coef(ij, 1)
         enddo
      enddo
      wf%fock_diagonal_cc2_ccs(:,1) = wf%fock_diagonal(:,1)
!
!     ::::::::::::::::::::::::::::::::
!     -::-   CCSD/CC2 partition   -::-
!     ::::::::::::::::::::::::::::::::
!
!     :: Construct canonical occupied and vacant density matrices ::    
!
      call allocator(density_o, wf%n_ao, wf%n_ao)
      call allocator(density_v, wf%n_ao, wf%n_ao)
      density_o = zero
      density_v = zero
!
      call wf%construct_density_matrices(density_o, density_v, wf%mo_coef, wf%n_o, wf%n_v)
!
!     No CCS region
!
      wf%n_CCS_o = 0
      wf%n_CCS_v = 0
!
      offset_o = 1
      offset_v = 1 + wf%n_o
!
!  
!     Get CC2/CCSD-active atoms
!
      n_CCSD_atoms =  get_number_of_active_atoms(unit_cholesky_decomp, 'CCSD ') 
!
!     Allocate active atoms list
!  
      call allocator_int(active_atoms, n_CCSD_atoms, 1)
!  
!     Get list of active atoms
!  
      call get_active_atoms(unit_cholesky_decomp, active_atoms,  n_CCSD_atoms, 'CCSD  ')
!  
!     Sanity check on active atoms
!  
      if (n_CCSD_atoms .gt. n_nuclei) then
         write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
         stop
      endif
!  
      do i = 1, n_CCSD_atoms
         if (active_atoms(i,1) .gt. n_nuclei) then
            write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
            stop
         endif
      enddo 
!  
!     :: Constructing active (CCSD) localized Cholesky orbitals ::
!  
      n_active_aos = 0
      do i = 1, n_CCSD_atoms
         n_active_aos = n_active_aos + n_ao_on_center(active_atoms(i,1),1)
      enddo
!  
!     Construct active_ao_index_list
!  
      call allocator_int(active_ao_index_list, n_active_aos, 1)
      call construct_active_ao_index_list(active_ao_index_list, n_active_aos, active_atoms, &
                                                  n_CCSD_atoms, ao_center_info, wf%n_ao) 
!  
!     Occupied part
!  
      n_vectors_o = 0 
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                              density_o, n_vectors_o,&
                                             .true., n_active_aos, active_ao_index_list)
!  
!     Virtual part
!  
      n_vectors_v = 0 
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                              density_v, n_vectors_v,&
                                              .true., n_active_aos, active_ao_index_list)
!  
!     Save active space information
!  
      wf%n_CCSD_o = n_vectors_o
      wf%n_CCSD_v = n_vectors_v
!
      wf%first_CCSD_o = offset_o
      wf%first_CCSD_v = offset_v - wf%n_o
!  
!     Calculate new offset         
!  
      offset_o = offset_o + n_vectors_o
      offset_v = offset_v + n_vectors_v
!  
      call deallocator_int(active_atoms, n_CCSD_atoms, 1)
      call deallocator_int(active_ao_index_list, n_active_aos, 1)
!
!     :: CC2  localized Cholesky orbitals  ::
!
!     Occupied part
!
      n_vectors_o = 0 
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                           density_o, n_vectors_o,&
                                           .false., n_active_aos)
!
!     Virtual part
!
      n_vectors_v = 0 
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                           density_v, n_vectors_v,&
                                           .false., n_active_aos)
!
!     Save inactive space information
!
      wf%n_CC2_o = n_vectors_o
      wf%n_CC2_v = n_vectors_v
!
      wf%first_CC2_o = offset_o
      wf%first_CC2_v = offset_v - wf%n_o
!
      wf%mo_coef       = orbitals
      wf%fock_diagonal = orbital_energies
!
!     Deallocations
!
      call deallocator(density_v, wf%n_ao, wf%n_ao)   
      call deallocator(density_o, wf%n_ao, wf%n_ao)
!
      call deallocator(orbitals, wf%n_ao, wf%n_mo)
      call deallocator(orbital_energies, wf%n_mo, 1) 
!
   end subroutine cholesky_localization_CCSD_CC2_mlccsd
!
!
   module subroutine construct_MO_transformation_matrix_mlccsd(wf)
!!
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(:,:), allocatable :: X
      real(dp), dimension(:,:), allocatable :: X2
      real(dp), dimension(:,:), allocatable :: S_packed, S
!
      integer(i15) :: a = 0, b = 0, i = 0, j = 0
!
      integer(i15) :: unit_overlap = -1, ioerror = 0
!
!
!     Read overlap matrix
!
      call generate_unit_identifier(unit_overlap)
      open(unit=unit_overlap, file='MLCC_OVERLAP', status='unknown', form='formatted', iostat=ioerror)
      if (ioerror .ne. 0) write(unit_output,*)'WARNING: Error while opening MLCC_OVERLAP'
      rewind(unit_overlap)
!
      call allocator(S_packed, wf%n_ao*(wf%n_ao+1)/2,1)
      read(unit_overlap,*)S_packed
      close(unit_overlap)
!
      call allocator(S, wf%n_ao, wf%n_ao)
!
      call squareup(S_packed, S, wf%n_ao)
!
      call deallocator(S_packed, wf%n_ao*(wf%n_ao+1)/2,1)
!
      call allocator(X, wf%n_mo, wf%n_ao)
!
      call dgemm('T', 'N',    &
                  wf%n_mo,    &
                  wf%n_ao,    &
                  wf%n_ao,    &
                  one,        &
                  wf%mo_coef, &
                  wf%n_ao,    &
                  S,          &
                  wf%n_ao,    &
                  zero,       &
                  X,          &
                  wf%n_mo)
!
      call deallocator(S, wf%n_ao, wf%n_ao)
      call allocator(X2, wf%n_mo, wf%n_mo)
!
      call dgemm('N', 'N',            &
                  wf%n_mo,            &
                  wf%n_mo,            &
                  wf%n_ao,            &
                  one,                &
                  X,                  &
                  wf%n_mo,            &
                  wf%mo_coef_cc2_ccs, &
                  wf%n_ao,            &
                  zero,               &
                  X2,                 &
                  wf%n_mo)
!
      call deallocator(X, wf%n_mo, wf%n_ao)
!
      call deallocator(wf%mo_coef_cc2_ccs, wf%n_ao, wf%n_mo)
      call allocator(wf%T_o, wf%n_o, wf%n_o)
      call allocator(wf%T_v, wf%n_v, wf%n_v)
!
!     Construct active occupied and active virtual transformation matrices
!
      do i = 1, wf%n_o
         do j = 1, wf%n_o
            wf%T_o(i,j) = X2(i, j)
         enddo
      enddo
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
            wf%T_v(a,b) = X2(wf%n_o + a, wf%n_o + b)
         enddo
      enddo
!
      call deallocator(X2, wf%n_mo, wf%n_mo)
!
   end subroutine construct_MO_transformation_matrix_mlccsd
!
!
   module subroutine cnto_orbital_drv_mlccsd(wf)
!!
!!    CNTO orbital driver,
!!    Written by Sarai D. Folkestad, June 2017.
!!
      implicit none 
!
      class(mlccsd) :: wf
!
      real(dp) :: start_cnto = 0, end_cnto = 0
!     
!     Prints
!     
      write(unit_output,'(/t3,a/)') ':: CNTO orbital partitioning for MLCCSD calculation '
!
!     Timings 
!
      call cpu_time(start_cnto)
!
!     CNTO orbital selection for CCSD
!
      if (wf%tasks%excited_state) then
!
         call wf%ccsd_cnto
!
      elseif (wf%tasks%core_excited_state) then
!
         call wf%ccsd_cnto_cvs
!
      else
!
         write(unit_output,*)'WARNING: cntos without excited state calculation makes no sense.'
      endif
      
!
!     Timings
!
      call cpu_time(end_cnto)
!
      if (timings) write(unit_output,'(/t3,a50,f14.8/)') 'Total time in cnto orbital construction (seconds):', end_cnto - start_cnto
      flush(unit_output)
!
   end subroutine cnto_orbital_drv_mlccsd
!
  module subroutine ccsd_cnto_mlccsd(wf)
!!
!!    CNTO orbital driver,
!!    Written by Sarai D. Folkestad, June 2017.
!!
!!    A CCS calculation ground state and excited states is performed.
!!    The M and N matrices are then constructed, 
!! 
!!       M_ij = sum_a R1_ai*R1_aj + sum_a R2_ai*R2_aj + ...
!!       N_ab = sum_i R1_ai*R1_bi + sum_a R2_ai*R2_bi + ...
!!   
!!    where Ri_ai is the i'th single excitation vector obtained from the CCS calculation. 
!!    The transformation matrices for the occupied and virtual part
!!    are constructed by diagonalizing M and N. The number of active occupied
!!    and virtual orbitals are determined from δ_o and δ_v
!!
!!       1 - sum_i λ^o_i < δ_o
!!       1 - sum_i λ^v_i < δ_v
!!
!!    Where the orbitals of highest eigenvalues λ^o/λ^v are selected first.
!!
!!    Fock matrix is block diagonalized in active and inactive blocks in order to obtain 
!!    the orbitals and orbital energies used in the CC2 calculation.
!!

      implicit none 
!
      class(mlccsd) :: wf
!
      type(mlcc2), allocatable :: cc2_wf
!
      integer(i15) :: unit_solution = -1
      integer(i15) :: ioerror = 0
      integer(i15) :: lower_level_n_singlet_states
      integer(i15) :: state
      integer(i15) :: cc2_n_x2am, cc2_n_parameters
!
      real(dp), dimension(:,:), allocatable :: solution
      real(dp), dimension(:,:), allocatable :: R
      real(dp), dimension(:,:), allocatable :: R_sum
      real(dp), dimension(:,:), allocatable :: R_sum_restricted
      real(dp), dimension(:,:), allocatable :: R_a_i
      real(dp), dimension(:,:), allocatable :: R_ai_bj
      real(dp), dimension(:,:), allocatable :: R_aib_k
      real(dp), dimension(:,:), allocatable :: R_a_icj
      real(dp), dimension(:,:), allocatable :: M_i_j
      real(dp), dimension(:,:), allocatable :: N_a_b
      real(dp), dimension(:,:), allocatable :: M
      real(dp), dimension(:,:), allocatable :: N
      real(dp), dimension(:,:), allocatable :: X
      real(dp), dimension(:,:), allocatable :: eigenvalues_o
      real(dp), dimension(:,:), allocatable :: eigenvalues_v
      real(dp), dimension(:,:), allocatable :: work
      real(dp), dimension(:,:), allocatable :: C_o
      real(dp), dimension(:,:), allocatable :: C_v
      real(dp), dimension(:,:), allocatable :: C_o_transformed
      real(dp), dimension(:,:), allocatable :: C_v_transformed
      real(dp), dimension(:,:), allocatable :: fock_o
      real(dp), dimension(:,:), allocatable :: fock_v
      real(dp), dimension(:,:), allocatable :: orbital_energies
      real(dp), dimension(:,:), allocatable :: test1
      real(dp), dimension(:,:), allocatable :: test2
!
      integer(i15) :: unit_dt = -1 , unit_t_dt = -1
!
      integer(i15) :: info
      integer(i15) :: i = 0, j = 0, a = 0, b = 0, c = 0, k = 0
      integer(i15) :: ai = 0, ij = 0, cj = 0, bj = 0, bk = 0, ai_full
      integer(i15) :: aibj = 0, aib = 0, icj = 0
!
      real(dp) :: trace, ddot, sum_o, sum_v, norm
      integer(i15) :: n_CC2_o, n_CC2_v
!
!     ::::::::::::::::::::::::::::::::::::::::::::::::
!     -::- Running lower level method calculation -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::::
!
!     Allocate lower level method
!
      allocate(cc2_wf)
!
!     Set calculation tasks
!
      cc2_wf%tasks%ground_state = .true.
      cc2_wf%tasks%excited_state = .true.
!
!     Set calculation settings
!
      cc2_wf%settings = wf%settings
! 
!     Set number of excitations to use for cnto generation
!
      cc2_wf%tasks%n_singlet_states = wf%tasks%n_singlet_states
!
!     Set convergence threshold for lower lying method
!
      cc2_wf%settings%energy_threshold = 1.0D-04 
      cc2_wf%settings%equation_threshold = 1.0D-04 
!
!     Set mlcc settings
!
!
      if(wf%mlcc_settings%CCS .and. wf%mlcc_settings%CC2) then
!
         cc2_wf%mlcc_settings%CCS      = .true.
         cc2_wf%mlcc_settings%cnto     = .true.
         cc2_wf%mlcc_settings%cholesky = .false.
!
      else
!
         cc2_wf%mlcc_settings%CC2      = .true.
         cc2_wf%mlcc_settings%CCS = .false.
!
      endif
!
!     Initialize lower level method
!  
      call cc2_wf%init
!
!     Call driver of lower level method
!
      call cc2_wf%drv
!
      cc2_n_parameters = cc2_wf%n_parameters
      cc2_n_x2am = cc2_wf%n_x2am
!
      do i = 1, wf%n_ao
         do j = 1, wf%n_mo
            ij = index_two(i, j, wf%n_ao)
            wf%mo_coef_cc2_ccs(i, j) = cc2_wf%mo_coef(ij, 1)
         enddo
      enddo
!
      wf%fock_diagonal_cc2_ccs = cc2_wf%fock_diagonal
!
      wf%mo_coef = cc2_wf%mo_coef
      wf%fock_diagonal = cc2_wf%fock_diagonal
!
      n_CC2_o = cc2_wf%n_CC2_o
      n_CC2_v = cc2_wf%n_CC2_v
!
!     Save info on CCS space
!
      wf%n_CCS_o = cc2_wf%n_CCS_o
      wf%n_CCS_v = cc2_wf%n_CCS_v
!
      wf%first_CCS_o = cc2_wf%first_CCS_o
      wf%first_CCS_v = cc2_wf%first_CCS_v
!
!     Deallocate lower level method
!
      deallocate(cc2_wf)
!
!     ::::::::::::::::::::::::::::::::::::::::::::::
!     -::- Construct CNTO transformation matrix -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::
!
      wf%mlcc_settings%delta_o = 1.0d-3
      wf%mlcc_settings%delta_v = 1.0d-5 ! THESE SHOULD BE SET IN INPUT
!
!     Open file of CC2 solution vectors
!
      call generate_unit_identifier(unit_solution)
!
      open(unit=unit_solution, file=wf%excited_state_task, action='read', status='unknown', &
        access='direct', form='unformatted', recl=dp*(cc2_n_parameters), iostat=ioerror) 
!
      if (ioerror .ne. 0) write(unit_output,*) 'Error while opening solution file', ioerror
!
      call allocator(R, cc2_n_parameters, 1)
      call allocator(R_sum, cc2_n_parameters, 1)
      R_sum = zero
!
!     Construct sum of excitation vectors
!
      do state = 1, wf%tasks%n_singlet_states
!
         read(unit=unit_solution, rec=state) R
!
         call daxpy(cc2_n_parameters, one, R, 1, R_sum, 1)
!
      enddo
!
!     Done with file, delete it
!
      close(unit_solution, status='delete')
!
      call deallocator(R, cc2_n_parameters, 1)
!
!
!     Restict R_sum such that only CC2 orbitals are included
!
      call allocator(R_sum_restricted, n_CC2_o*n_CC2_v + cc2_n_x2am, 1)
      do a = 1, n_CC2_v
         do i = 1, n_CC2_o
            ai = index_two(a, i, n_CC2_v)
            ai_full = index_two(a, i, wf%n_v)
            R_sum_restricted(ai, 1) = R_sum(ai_full, 1)
         enddo
      enddo
!
      do a = 1, n_CC2_v
         do i = 1, n_CC2_o
            ai = index_two(a, i, n_CC2_v)
            do b = 1, n_CC2_v
               do j = 1, n_CC2_o
                  bj = index_two(b, j, n_CC2_v)
                  aibj = index_packed(ai, bj)
                  R_sum_restricted((n_CC2_o)*(n_CC2_v) + aibj, 1) &
                        = R_sum((wf%n_o)*(wf%n_v) + aibj, 1)
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(R_sum, cc2_n_parameters, 1)
!
      norm = sqrt(ddot(cc2_n_x2am + n_CC2_o*n_CC2_v, R_sum_restricted, 1, R_sum_restricted, 1))
      call dscal(cc2_n_x2am + n_CC2_o*n_CC2_v, one/norm, R_sum_restricted, 1)
!
      call allocator(R_a_i, n_cc2_v, n_cc2_o)
!
      do a = 1, n_cc2_v
         do i = 1, n_cc2_o
            ai = index_two(a, i, n_CC2_v)
            R_a_i(a, i) = R_sum_restricted(ai, 1)
         enddo
      enddo
!
      call allocator(R_ai_bj, n_cc2_o*n_cc2_v, n_cc2_o*n_cc2_v)
!
      do ai = 1, n_cc2_o*n_cc2_v
         do bj = 1, n_cc2_o*n_cc2_v
            aibj = index_packed(ai, bj)
            R_ai_bj(ai, bj) = R_sum_restricted(aibj + n_cc2_o*n_cc2_v, 1)
         enddo
      enddo
!
      call deallocator(R_sum_restricted, n_CC2_o*n_CC2_v + cc2_n_x2am, 1)

!
!     Construct M and N
!
      call allocator(M_i_j, n_CC2_o, n_CC2_o)
      call allocator(N_a_b, n_CC2_v, n_CC2_v)
!
!     Singles contribution
!
      call dgemm('T', 'N',    &
                     n_CC2_o, &
                     n_CC2_o, &
                     n_CC2_v, &
                     one,     &
                     R_a_i,   &
                     n_CC2_v, &
                     R_a_i,   &
                     n_CC2_v, &
                     zero,    &
                     M_i_j,   &
                     n_CC2_o)
!
      call dgemm('N', 'T', &
                     n_CC2_v, &
                     n_CC2_v, &
                     n_CC2_o, &
                     one,     &
                     R_a_i,   &
                     n_CC2_v, &
                     R_a_i,   &
                     n_CC2_v, &
                     zero,    &
                     N_a_b,   &
                     n_CC2_v)
!
      call deallocator(R_a_i, n_cc2_v, n_cc2_o)
!
      call allocator(R_aib_k, (n_CC2_v**2)*(n_CC2_o), n_CC2_o)
!
      do a = 1, n_CC2_v
         do b = 1, n_CC2_v
            do k = 1, n_CC2_o
               do i = 1, n_CC2_o
                  aib = index_three(a, i, b, n_CC2_v, n_CC2_o)
                  ai = index_two(a, i, n_CC2_v)
                  bk = index_two(b, k, n_CC2_v)
                  R_aib_k(aib, k) = R_ai_bj(ai, bk)
               enddo
            enddo
         enddo
      enddo
!
      call dgemm('T', 'N',                   &
                     n_CC2_o,                &
                     n_CC2_o,                &
                     (n_CC2_v**2)*(n_CC2_o), &
                     half,                   &
                     R_aib_k,                & ! R_ai,bk
                     (n_CC2_v**2)*(n_CC2_o), &
                     R_aib_k,                & ! R_aj,bk
                     (n_CC2_v**2)*(n_CC2_o), &
                     one,                    &
                     M_i_j,                  &
                     n_CC2_o)
!
      call deallocator(R_aib_k, (n_CC2_v**2)*(n_CC2_o), n_CC2_o)
!
      do i = 1, n_CC2_o
         do a = 1, n_CC2_v
            ai = index_two(a, i, n_CC2_v)
            M_i_j(i, i) = M_i_j(i, i) + half*R_ai_bj(ai, ai)
         enddo
      enddo
!
      call allocator(R_a_icj, n_CC2_v, (n_CC2_o**2)*(n_CC2_v))
!
      do a = 1, n_CC2_v
         do c = 1, n_CC2_v
            do j = 1, n_CC2_o
               do i = 1, n_CC2_o
                  icj = index_three(i, c, j, n_CC2_o, n_CC2_v)
                  ai = index_two(a, i, n_CC2_v)
                  cj = index_two(c, j, n_CC2_v)
                  R_a_icj(a, icj) = R_ai_bj(ai, cj)
               enddo
            enddo
         enddo
      enddo
      call dgemm('N', 'T',                   &
                     n_CC2_v,                &
                     n_CC2_v,                &
                     (n_CC2_o**2)*(n_CC2_v), &
                     half,                   &
                     R_a_icj,                & ! R_ai,cj
                     n_CC2_v,                &
                     R_a_icj,                & ! R_bi,cj
                     n_CC2_v,                &
                     one,                    &
                     N_a_b,                  &
                     n_CC2_v)
!
      call deallocator(R_a_icj, n_CC2_v, (n_CC2_o**2)*(n_CC2_v))
!
      do a = 1, n_CC2_v
         do i = 1, n_CC2_o
            ai = index_two(a, i, n_CC2_v)
            N_a_b(a, a) = N_a_b(a, a) + half*R_ai_bj(ai, ai)
         enddo
      enddo
!
      call deallocator(R_ai_bj, n_cc2_o*n_cc2_v, n_cc2_o*n_cc2_v)
!
!
!     :: Diagonalize M and N matrix ::
!
      call allocator(eigenvalues_o, n_cc2_o, 1)
      call allocator(work, 4*(n_cc2_o), 1)
      work = zero
!
      call dsyev('V','U', &
                  n_cc2_o, &
                  M_i_j, &
                  n_cc2_o, &
                  eigenvalues_o, &
                  work, & 
                  4*(n_cc2_o), &
                  info)
!
      call deallocator(work, 4*(n_cc2_o), 1)
!
      call allocator(eigenvalues_v, n_cc2_v, 1)
      call allocator(work, 4*(n_cc2_v), 1)
      work = zero
!
      call dsyev('V','U', &
                  n_cc2_v, &
                  N_a_b, &
                  n_cc2_v, &
                  eigenvalues_v, &
                  work, & 
                  4*(n_cc2_v), &
                  info)
!
      call deallocator(work, 4*(n_cc2_v), 1)
!
!     :: Reorder M and N ::
!
!     dsyev orderes eigenvalues and corresponding eigenvectors in ascending order.
!     We wish to select active space according to highest eigenvalues, thus we must reorder
!
      call allocator(M, n_cc2_o, n_cc2_o)
!
      do i = 1, n_cc2_o
!
         j = i -1
         M(:,i) = M_i_j(:, n_cc2_o - j)
!
      enddo
!
      call deallocator(M_i_j, n_CC2_o, n_CC2_o)
!
      call allocator(N, n_cc2_v, n_cc2_v)
!
      do a = 1, n_cc2_v
!
         b = a -1
         N(:,a) = N_a_b(:,n_cc2_v - b)
!
      enddo
!
      call deallocator(N_a_b, n_CC2_v, n_CC2_v)
!
!
!     Transform C to CNTO
!
      call allocator(C_o, wf%n_ao, n_cc2_o)
      call allocator(C_v, wf%n_ao, n_cc2_v)
      C_o = zero
      C_v = zero      
!
     do i = 1, wf%n_ao
!
        do j = 1, n_cc2_o
           ij = index_two(i, j, wf%n_ao)
           C_o(i, j) = wf%mo_coef(ij, 1) 
        enddo
!
        do j = 1, n_cc2_v 
           ij = index_two(i, j + wf%n_o, wf%n_ao)
           C_v(i, j) = wf%mo_coef(ij, 1) 
        enddo
!
     enddo
!
      call allocator(C_o_transformed, wf%n_ao, n_cc2_o)
      call dgemm('N', 'N',    &
                  wf%n_ao,    &
                  n_cc2_o,    &
                  n_cc2_o,    &
                  one,        &
                  C_o,        &
                  wf%n_ao,    &
                  M,          &
                  n_cc2_o,    &
                  zero,       &
                  C_o_transformed, &
                  wf%n_ao)
!
      call deallocator(C_o, wf%n_ao, n_cc2_o)
!
      call allocator(C_v_transformed, wf%n_ao, n_cc2_v)
      call dgemm('N', 'N',                   &
                  wf%n_ao,                   &
                  n_cc2_v,                   &
                  n_cc2_v,                   &
                  one,                       &
                  C_v,                       &
                  wf%n_ao,                   &
                  N,                         &
                  n_cc2_v,                   &
                  zero,                      &
                  C_v_transformed,           &
                  wf%n_ao)
!
      call deallocator(C_v, wf%n_ao, n_cc2_v)
      call deallocator(N, n_cc2_v, n_cc2_v)
      call deallocator(M, n_cc2_o, n_cc2_o)  
!
      do i = 1, wf%n_ao
!

         do j = 1, n_cc2_o
            ij = index_two(i, j, wf%n_ao)
            wf%mo_coef(ij, 1) = C_o_transformed(i, j) 
         enddo
!
         do j = 1, n_cc2_v 
            ij = index_two(i, j + wf%n_o, wf%n_ao)
            wf%mo_coef(ij, 1) = C_v_transformed(i, j) 
         enddo
!
      enddo
!
      call deallocator(C_o_transformed, wf%n_ao, n_cc2_o)
      call deallocator(C_v_transformed, wf%n_ao, n_cc2_v)
!
!     :: Determine number of active orbitals ::
!
      sum_o       = 1
      wf%n_CCSD_o = 1
!
      do while ((sum_o .gt. wf%mlcc_settings%delta_o) .and. (wf%n_CCSD_o .lt. n_cc2_o))
!
         sum_o = sum_o - eigenvalues_o(n_CC2_o - (wf%n_CCSD_o - 1), 1)
         wf%n_CCSD_o = wf%n_CCSD_o + 1
!
      enddo

!
      sum_v      = 1
      wf%n_CCSD_v = 1
!
      do while (sum_v .gt. wf%mlcc_settings%delta_v .and. (wf%n_CCSD_v .lt. n_cc2_v))
!
         sum_v = sum_v - eigenvalues_v(n_CC2_v - (wf%n_CCSD_v - 1), 1)
         wf%n_CCSD_v = wf%n_CCSD_v + 1
!
      enddo 
!
!     Save information to object
!
      wf%first_CCSD_o = 1
      wf%first_CCSD_v = 1
!
      wf%n_CC2_o = n_cc2_o - wf%n_CCSD_o
      wf%n_CC2_v = n_cc2_v - wf%n_CCSD_v
!
      wf%first_CC2_o = 1 + wf%n_CCSD_o
      wf%first_CC2_v = 1 + wf%n_CCSD_v
!
      call deallocator(eigenvalues_o, n_cc2_o, 1)
      call deallocator(eigenvalues_v, n_cc2_v, 1)
      call wf%print_cnto_info
!     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     -::- Finding orbital energies and new block diagonal C matrix -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     Initialize amplitudes and associated attributes
!
      call wf%initialize_amplitudes
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf%read_transform_cholesky
!
!     Initialize fock matrix
!
      call wf%initialize_fock_matrix
!
!     :: Occupied Orbitals ::
!
!     Diagonalize active-active block
!
      if (wf%n_CCSD_o .gt. 0) then
!
         call allocator(work, 4*(wf%n_CCSD_o), 1)
         call allocator(orbital_energies, (wf%n_CCSD_o), 1)
         work = zero
!
         call dsyev('V','U',              &
                     (wf%n_CCSD_o),       &
                     wf%fock_ij,          &
                     wf%n_o,              &
                     orbital_energies,    &
                     work,                & 
                     4*(wf%n_CCSD_o),     &
                     info)
!
         call deallocator(work, 4*(wf%n_CCSD_o), 1)
!
         if (info .ne. 0) then
            write(unit_output,*)'WARNING: Diagonalization of active virtual block not successful. '
            stop
         endif
!
!        Setting orbital energies
!
         do j = 1, wf%n_CCSD_o
!
            wf%fock_diagonal(j,1) = orbital_energies(j,1)
!
         enddo
!
         call deallocator(orbital_energies, (wf%n_CCSD_o), 1)
      endif
!
!     Diagonalize inactive-inactive block 
!
      if (wf%n_CC2_o .gt. 0) then
         call allocator(work, 4*wf%n_CC2_o, 1)
         call allocator(orbital_energies, wf%n_CC2_o, 1)
         orbital_energies = zero
         work = zero
!
         call dsyev('V','U',                                   &
                     wf%n_CC2_o,                               &
                     wf%fock_ij(wf%first_CC2_o, wf%first_CC2_o),&
                     wf%n_o,                                   &
                     orbital_energies,                         &
                     work,                                     & 
                     4*wf%n_CC2_o,                             &
                     info)
!
         call deallocator(work, 4*wf%n_CC2_o, 1)
!
         if (info .ne. 0) then
            write(unit_output,*)'WARNING: Diagonalization of inactive virtual block not successful.'
            stop
         endif
!
!        Setting orbital energies
!
         do j = 1, wf%n_CC2_o
!
            wf%fock_diagonal(j + wf%first_CC2_o - 1,1) = orbital_energies(j,1)
!
         enddo
         call deallocator(orbital_energies, wf%n_CC2_o, 1)
      endif
!
!     Transform C-matrix (active occupied block)
!
      call allocator(C_o, wf%n_ao, wf%n_o)
      C_o = zero     
!  
      do i = 1, wf%n_ao
!  
         do j = 1, wf%n_o
           ij = index_two(i, j, wf%n_ao)
           C_o(i, j) = wf%mo_coef(ij, 1) 
         enddo
!  
      enddo
!
      if (wf%n_CCSD_o .gt. 0) then
!  
         call allocator(C_o_transformed, wf%n_ao, wf%n_CCSD_o)
!  
         call dgemm('N', 'N',          &
                     wf%n_ao,          &
                     wf%n_CCSD_o,      &
                     wf%n_CCSD_o,      &
                     one,              &
                     C_o,              &
                     wf%n_ao,          &
                     wf%fock_ij,       &
                     wf%n_o,           &
                     zero,             &
                     C_o_transformed,  &
                     wf%n_ao)
!  
         do i = 1, wf%n_ao
!  
            do j = 1, wf%n_CCSD_o
               ij = index_two(i, j, wf%n_ao)
               wf%mo_coef(ij, 1) = C_o_transformed(i, j) 
            enddo
!  
         enddo
!  
         call deallocator(C_o_transformed, wf%n_ao, wf%n_CCSD_o)
      endif
!
!     Transform C-matrix (inactive occupied block)
!
      if (wf%n_CC2_o .gt. 0) then
         call allocator(C_o_transformed, wf%n_ao, wf%n_CC2_o)
         call dgemm('N', 'N',    &
                     wf%n_ao,    &
                     wf%n_CC2_o, &
                     wf%n_CC2_o, &
                     one,        &
                     C_o(1, wf%first_CC2_o),        &
                     wf%n_ao,    &
                     wf%fock_ij(wf%first_CC2_o, wf%first_CC2_o), &
                     wf%n_o,     &
                     zero,       &
                     C_o_transformed, &
                     wf%n_ao)
!
         call deallocator(C_o, wf%n_ao, wf%n_o)
!
         do i = 1, wf%n_ao
!
            do j = 1, wf%n_CC2_o
               ij = index_two(i, j + wf%first_CC2_o -1, wf%n_ao)
               wf%mo_coef(ij, 1) = C_o_transformed(i, j) 
            enddo
!
         enddo
!
         call deallocator(C_o_transformed, wf%n_ao, wf%n_CC2_o)
      endif
!
!     :: Vacant orbitals ::
!
!     Diagonalize active-active block
!
      if (wf%n_CCSD_v .gt. 0) then
         call allocator(work, 4*(wf%n_CCSD_v), 1)
         call allocator(orbital_energies, (wf%n_CCSD_v), 1)
         work = zero
!
         call dsyev('V','U',              &
                     (wf%n_CCSD_v),   &
                     wf%fock_ab,          &
                     wf%n_v,              &
                     orbital_energies,    &
                     work,                & 
                     4*(wf%n_CCSD_v), &
                     info)
!
         call deallocator(work, 4*(wf%n_CCSD_v), 1)
!
         if (info .ne. 0) then
            write(unit_output,*)'WARNING: Diagonalization of active virtual block not successful. '
                  stop
         endif
!
!        Setting orbital energies
!
         do j = 1, wf%n_CCSD_v
!
            wf%fock_diagonal(j + wf%n_o ,1) = orbital_energies(j,1)
!
         enddo
!
         call deallocator(orbital_energies, (wf%n_CCSD_v), 1)
!
      endif
!
!     Diagonalize inactive-inactive block 
!
      if (wf%n_CC2_v .gt. 0) then
         call allocator(work, 4*wf%n_CC2_v, 1)
         call allocator(orbital_energies, wf%n_CC2_v, 1)
         orbital_energies = zero
         work = zero
!
         call dsyev('V','U',                                      &
                     wf%n_CC2_v,                                  &
                     wf%fock_ab(wf%first_CC2_v, wf%first_CC2_v),  &
                     wf%n_v,                                      &
                     orbital_energies,                            &
                     work,                                        & 
                     4*(wf%n_CC2_v),                              &
                     info)
!
         call deallocator(work, 4*wf%n_CC2_v, 1)
!
         if (info .ne. 0) then
            write(unit_output,*)'WARNING: Diagonalization of inactive virtual block not successful.'
            stop
         endif
!
!        Setting orbital energies
!
         do j = 1, wf%n_CC2_v
!
            wf%fock_diagonal(j + wf%n_o + wf%first_CC2_v - 1,1) = orbital_energies(j,1)
!
         enddo

         call deallocator(orbital_energies, wf%n_CC2_v, 1)
!
      endif
!
!     Transform C-matrix (active virtual block)
!
      call allocator(C_v, wf%n_ao, wf%n_v)
      C_v = zero     
!
      do i = 1, wf%n_ao
!
         do j = 1, wf%n_v
            ij = index_two(i, j + wf%n_o, wf%n_ao)
            C_v(i, j) = wf%mo_coef(ij, 1) 
         enddo
!
      enddo
!
      if (wf%n_CCSD_v .gt. 0) then
!
         call allocator(C_v_transformed, wf%n_ao, wf%n_CCSD_v)
         call dgemm('N', 'N',    &
                     wf%n_ao,    &
                     wf%n_CCSD_v,&
                     wf%n_CCSD_v,&
                     one,        &
                     C_v,        &
                     wf%n_ao,    &
                     wf%fock_ab, &
                     wf%n_v,     &
                     zero,       &
                     C_v_transformed, &
                     wf%n_ao)
!
         do i = 1, wf%n_ao
            do j = 1, wf%n_CCSD_v
               ij = index_two(i, j + wf%n_o, wf%n_ao)
               wf%mo_coef(ij, 1) = C_v_transformed(i, j) 
            enddo
         enddo
!
         call deallocator(C_v_transformed, wf%n_ao, wf%n_CCSD_v)
!
      endif
!
!     Transform C-matrix (inactive virtual block)
!
      if (wf%n_CC2_v .gt. 0) then
!
         call allocator(C_v_transformed, wf%n_ao, wf%n_CC2_v)
         call dgemm('N', 'N',    &
                     wf%n_ao,    &
                     wf%n_CC2_v, &
                     wf%n_CC2_v, &
                     one,        &
                     C_v(1, wf%first_CC2_v),        &
                     wf%n_ao,    &
                     wf%fock_ab(wf%first_CC2_v, wf%first_CC2_v), &
                     wf%n_v,     &
                     zero,       &
                     C_v_transformed, &
                     wf%n_ao)
!
         do i = 1, wf%n_ao

            do j = 1, wf%n_CC2_v
               ij = index_two(i, j + wf%n_o + wf%first_CC2_v - 1, wf%n_ao)
               wf%mo_coef(ij, 1) = C_v_transformed(i, j) 
            enddo
         enddo
!
         call deallocator(C_v_transformed, wf%n_ao, wf%n_CC2_v)
!
      endif
!
      call deallocator(C_v, wf%n_ao, wf%n_v)
      
!
   end subroutine ccsd_cnto_mlccsd
!
!
 module subroutine ccsd_cnto_cvs_mlccsd(wf)
!!
!!    CNTO orbital driver,
!!    Written by Sarai D. Folkestad, June 2017.
!!
!!    A CCS calculation ground state and excited states is performed.
!!    The M and N matrices are then constructed, 
!! 
!!       M_ij = sum_a R1_ai*R1_aj + sum_a R2_ai*R2_aj + ...
!!       N_ab = sum_i R1_ai*R1_bi + sum_a R2_ai*R2_bi + ...
!!   
!!    where Ri_ai is the i'th single excitation vector obtained from the CCS calculation. 
!!    The transformation matrices for the occupied and virtual part
!!    are constructed by diagonalizing M and N. The number of active occupied
!!    and virtual orbitals are determined from δ_o and δ_v
!!
!!       1 - sum_i λ^o_i < δ_o
!!       1 - sum_i λ^v_i < δ_v
!!
!!    Where the orbitals of highest eigenvalues λ^o/λ^v are selected first.
!!
!!    Fock matrix is block diagonalized in active and inactive blocks in order to obtain 
!!    the orbitals and orbital energies used in the CC2 calculation.
!!

      implicit none 
!
      class(mlccsd) :: wf
!
      type(mlcc2), allocatable :: cc2_wf
!
      integer(i15) :: unit_solution = -1
      integer(i15) :: ioerror = 0
      integer(i15) :: lower_level_n_singlet_states
      integer(i15) :: state
      integer(i15) :: cc2_n_x2am, cc2_n_parameters
!
      real(dp), dimension(:,:), allocatable :: solution
      real(dp), dimension(:,:), allocatable :: R
      real(dp), dimension(:,:), allocatable :: R_sum
      real(dp), dimension(:,:), allocatable :: R_sum_restricted
      real(dp), dimension(:,:), allocatable :: R_a_i
      real(dp), dimension(:,:), allocatable :: R_ai_bj
      real(dp), dimension(:,:), allocatable :: R_aib_k
      real(dp), dimension(:,:), allocatable :: R_a_icj
      real(dp), dimension(:,:), allocatable :: M_i_j
      real(dp), dimension(:,:), allocatable :: N_a_b
      real(dp), dimension(:,:), allocatable :: M
      real(dp), dimension(:,:), allocatable :: N
      real(dp), dimension(:,:), allocatable :: X
      real(dp), dimension(:,:), allocatable :: eigenvalues_o
      real(dp), dimension(:,:), allocatable :: eigenvalues_v
      real(dp), dimension(:,:), allocatable :: work
      real(dp), dimension(:,:), allocatable :: C_o
      real(dp), dimension(:,:), allocatable :: C_v
      real(dp), dimension(:,:), allocatable :: C_o_transformed
      real(dp), dimension(:,:), allocatable :: C_v_transformed
      real(dp), dimension(:,:), allocatable :: fock_o
      real(dp), dimension(:,:), allocatable :: fock_v
      real(dp), dimension(:,:), allocatable :: orbital_energies
      real(dp), dimension(:,:), allocatable :: test1
      real(dp), dimension(:,:), allocatable :: test2
!
      integer(i15) :: unit_dt = -1 , unit_t_dt = -1
!
      integer(i15) :: info
      integer(i15) :: i = 0, j = 0, a = 0, b = 0, c = 0, k = 0
      integer(i15) :: ai = 0, ij = 0, cj = 0, bj = 0, bk = 0, ai_full
      integer(i15) :: aibj = 0, aibk = 0, aib = 0, icj = 0
!
      real(dp) :: trace, ddot, sum_o, sum_v, norm
      integer(i15) :: n_CC2_o, n_CC2_v
!
!     ::::::::::::::::::::::::::::::::::::::::::::::::
!     -::- Running lower level method calculation -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::::
!
!     Allocate lower level method
!
      allocate(cc2_wf)
!
!     Set calculation tasks
!
      cc2_wf%tasks%ground_state = .true.
      cc2_wf%tasks%core_excited_state = .true.
      cc2_wf%tasks%n_cores = wf%tasks%n_cores
      call allocator_int(cc2_wf%tasks%cores, cc2_wf%tasks%n_cores, 1)
      cc2_wf%tasks%cores = wf%tasks%cores
!
!     Set calculation settings
!
      cc2_wf%settings = wf%settings
! 
!     Set number of excitations to use for cnto generation
!
      cc2_wf%tasks%n_singlet_states = wf%tasks%n_singlet_states
!
!     Set convergence threshold for lower lying method
!
      cc2_wf%settings%energy_threshold = 1.0D-04 
      cc2_wf%settings%equation_threshold = 1.0D-04 
!
!     Set mlcc settings
!
!
      if(wf%mlcc_settings%CCS .and. wf%mlcc_settings%CC2) then
!
         write(unit_output,*) 'Not implemented yet'
         stop
         cc2_wf%mlcc_settings%CCS      = .true.
         cc2_wf%mlcc_settings%cnto     = .true.
         cc2_wf%mlcc_settings%cholesky = .false.
!
      else
!
         cc2_wf%mlcc_settings%CC2      = .true.
         cc2_wf%mlcc_settings%CCS = .false.
!
      endif
!
!     Initialize lower level method
!  
      call cc2_wf%init
!
!     Call driver of lower level method
!
      call cc2_wf%drv
!
      cc2_n_parameters = cc2_wf%n_parameters
      cc2_n_x2am = cc2_wf%n_x2am
!
      do i = 1, wf%n_ao
         do j = 1, wf%n_mo
            ij = index_two(i, j, wf%n_ao)
            wf%mo_coef_cc2_ccs(i, j) = cc2_wf%mo_coef(ij, 1)
         enddo
      enddo
!
      wf%fock_diagonal_cc2_ccs = cc2_wf%fock_diagonal
!
      wf%mo_coef = cc2_wf%mo_coef
      wf%fock_diagonal = cc2_wf%fock_diagonal
!
      n_CC2_o = cc2_wf%n_CC2_o
      n_CC2_v = cc2_wf%n_CC2_v
!
!     Save info on CCS space
!
      wf%n_CCS_o = cc2_wf%n_CCS_o
      wf%n_CCS_v = cc2_wf%n_CCS_v
!
      wf%first_CCS_o = cc2_wf%first_CCS_o
      wf%first_CCS_v = cc2_wf%first_CCS_v
!
!     Deallocate lower level method
!
      deallocate(cc2_wf)
!
!     ::::::::::::::::::::::::::::::::::::::::::::::
!     -::- Construct CNTO transformation matrix -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::
!
      wf%mlcc_settings%delta_o = 1.0d-3
      wf%mlcc_settings%delta_v = 1.0d-5 ! THESE SHOULD BE SET IN INPUT
!
!     Open file of CC2 solution vectors
!
      call generate_unit_identifier(unit_solution)
!
      open(unit=unit_solution, file='right_core', action='read', status='unknown', &
        access='direct', form='unformatted', recl=dp*(cc2_n_parameters), iostat=ioerror) 
!
      if (ioerror .ne. 0) write(unit_output,*) 'Error while opening solution file', ioerror
!
      call allocator(R, cc2_n_parameters, 1)
      call allocator(R_sum, cc2_n_parameters, 1)
      R_sum = zero
!
!     Construct sum of excitation vectors
!
      do state = 1, wf%tasks%n_singlet_states
!
         read(unit=unit_solution, rec=state) R
!
         call daxpy(cc2_n_parameters, one, R, 1, R_sum, 1)
!
      enddo
!
!     Done with file, delete it
!
      close(unit_solution, status='delete')
!
      call deallocator(R, cc2_n_parameters, 1)
!
      call allocator(R_sum_restricted, n_CC2_o*n_CC2_v + cc2_n_x2am, 1)
      do a = 1, n_CC2_v
         do i = 1, n_CC2_o
            ai = index_two(a, i, n_CC2_v)
            ai_full = index_two(a, i, wf%n_v)
            R_sum_restricted(ai, 1) = R_sum(ai_full, 1)
         enddo
      enddo
!
      do a = 1, n_CC2_v
         do i = 1, n_CC2_o
            ai = index_two(a, i, n_CC2_v)
            do b = 1, n_CC2_v
               do j = 1, n_CC2_o
                  bj = index_two(b, j, n_CC2_v)
                  aibj = index_packed(ai, bj)
                  R_sum_restricted((n_CC2_o)*(n_CC2_v) + aibj, 1) &
                        = R_sum((wf%n_o)*(wf%n_v) + aibj, 1)
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(R_sum, cc2_n_parameters, 1)
!
      norm = sqrt(ddot(cc2_n_x2am + n_CC2_o*n_CC2_v, R_sum_restricted, 1, R_sum_restricted, 1))
      call dscal(cc2_n_x2am + n_CC2_o*n_CC2_v, one/norm, R_sum_restricted, 1)
!
      call allocator(R_a_i, n_cc2_v, n_cc2_o)
!
      do a = 1, n_cc2_v
         do i = 1, n_cc2_o
            ai = index_two(a, i, n_CC2_v)
            R_a_i(a, i) = R_sum_restricted(ai, 1)
         enddo
      enddo
!
      call allocator(R_ai_bj, n_cc2_o*n_cc2_v, n_cc2_o*n_cc2_v)
!
      do ai = 1, n_cc2_o*n_cc2_v
         do bj = 1, n_cc2_o*n_cc2_v
            aibj = index_packed(ai, bj)
            R_ai_bj(ai, bj) = R_sum_restricted(aibj + n_cc2_o*n_cc2_v, 1)
         enddo
      enddo
!
      call deallocator(R_sum_restricted, n_CC2_o*n_CC2_v + cc2_n_x2am, 1)

!
!     Construct M and N
!
      call allocator(M_i_j, n_CC2_o, n_CC2_o)
      call allocator(N_a_b, n_CC2_v, n_CC2_v)
!
!     Singles contribution
!
      call dgemm('T', 'N',    &
                     n_CC2_o, &
                     n_CC2_o, &
                     n_CC2_v, &
                     one,     &
                     R_a_i,   &
                     n_CC2_v, &
                     R_a_i,   &
                     n_CC2_v, &
                     zero,    &
                     M_i_j,   &
                     n_CC2_o)
!
      call dgemm('N', 'T', &
                     n_CC2_v, &
                     n_CC2_v, &
                     n_CC2_o, &
                     one,     &
                     R_a_i,   &
                     n_CC2_v, &
                     R_a_i,   &
                     n_CC2_v, &
                     zero,    &
                     N_a_b,   &
                     n_CC2_v)
!
!
      call deallocator(R_a_i, n_cc2_v, n_cc2_o)
!
      call allocator(R_aib_k, (n_CC2_v**2)*(n_CC2_o), n_CC2_o)
!
      do a = 1, n_CC2_v
         do b = 1, n_CC2_v
            do k = 1, n_CC2_o
               do i = 1, n_CC2_o
                  aib = index_three(a, i, b, n_CC2_v, n_CC2_o)
                  ai = index_two(a, i, n_CC2_v)
                  bk = index_two(b, k, n_CC2_v)
                  R_aib_k(aib, k) = R_ai_bj(ai, bk)
               enddo
            enddo
         enddo
      enddo
!
      call dgemm('T', 'N',                   &
                     n_CC2_o,                &
                     n_CC2_o,                &
                     (n_CC2_v**2)*(n_CC2_o), &
                     half,                   &
                     R_aib_k,                & ! R_ai,bk
                     (n_CC2_v**2)*(n_CC2_o), &
                     R_aib_k,                & ! R_aj,bk
                     (n_CC2_v**2)*(n_CC2_o), &
                     one,                    &
                     M_i_j,                  &
                     n_CC2_o)
!
      call deallocator(R_aib_k, (n_CC2_v**2)*(n_CC2_o), n_CC2_o)
!
      do i = 1, n_CC2_o
         do a = 1, n_CC2_v
            ai = index_two(a, i, n_CC2_v)
            M_i_j(i, i) = M_i_j(i, i) + half*R_ai_bj(ai, ai)
         enddo
      enddo
!
      call allocator(R_a_icj, n_CC2_v, (n_CC2_o**2)*(n_CC2_v))
!
      do a = 1, n_CC2_v
         do c = 1, n_CC2_v
            do j = 1, n_CC2_o
               do i = 1, n_CC2_o
                  icj = index_three(i, c, j, n_CC2_o, n_CC2_v)
                  ai = index_two(a, i, n_CC2_v)
                  cj = index_two(c, j, n_CC2_v)
                  R_a_icj(a, icj) = R_ai_bj(ai, cj)
               enddo
            enddo
         enddo
      enddo
      call dgemm('N', 'T',                   &
                     n_CC2_v,                &
                     n_CC2_v,                &
                     (n_CC2_o**2)*(n_CC2_v), &
                     half,                   &
                     R_a_icj,                & ! R_ai,cj
                     n_CC2_v,                &
                     R_a_icj,                & ! R_bi,cj
                     n_CC2_v,                &
                     one,                    &
                     N_a_b,                  &
                     n_CC2_v)
!
      call deallocator(R_a_icj, n_CC2_v, (n_CC2_o**2)*(n_CC2_v))
!
      do a = 1, n_CC2_v
         do i = 1, n_CC2_o
            ai = index_two(a, i, n_CC2_v)
            N_a_b(a, a) = N_a_b(a, a) + half*R_ai_bj(ai, ai)
         enddo
      enddo
!
      call deallocator(R_ai_bj, n_cc2_o*n_cc2_v, n_cc2_o*n_cc2_v)
!
!
!     :: Diagonalize M and N matrix ::
!
      call allocator(eigenvalues_o, n_cc2_o, 1)
      call allocator(work, 4*(n_cc2_o), 1)
      work = zero
!
      call dsyev('V','U', &
                  n_cc2_o, &
                  M_i_j, &
                  n_cc2_o, &
                  eigenvalues_o, &
                  work, & 
                  4*(n_cc2_o), &
                  info)
!
      call deallocator(work, 4*(n_cc2_o), 1)
!
      call allocator(eigenvalues_v, n_cc2_v, 1)
      call allocator(work, 4*(n_cc2_v), 1)
      work = zero
!
      call dsyev('V','U', &
                  n_cc2_v, &
                  N_a_b, &
                  n_cc2_v, &
                  eigenvalues_v, &
                  work, & 
                  4*(n_cc2_v), &
                  info)
!
      call deallocator(work, 4*(n_cc2_v), 1)
      write(unit_output,*)'9'
      flush(unit_output)
!
!     :: Reorder M and N ::
!
!     dsyev orderes eigenvalues and corresponding eigenvectors in ascending order.
!     We wish to select active space according to highest eigenvalues, thus we must reorder
!
      call allocator(M, n_cc2_o, n_cc2_o)
!
      do i = 1, n_cc2_o
!
         j = i -1
         M(:,i) = M_i_j(:, n_cc2_o - j)
!
      enddo
!
      call deallocator(M_i_j, n_CC2_o, n_CC2_o)
!
      call allocator(N, n_cc2_v, n_cc2_v)
!
      do a = 1, n_cc2_v
!
         b = a -1
         N(:,a) = N_a_b(:,n_cc2_v - b)
!
      enddo
!
      call deallocator(N_a_b, n_CC2_v, n_CC2_v)
!
!
!     Transform C to CNTO
!
      call allocator(C_o, wf%n_ao, n_cc2_o)
      call allocator(C_v, wf%n_ao, n_cc2_v)
      C_o = zero
      C_v = zero      
!
     do i = 1, wf%n_ao
!
        do j = 1, n_cc2_o
           ij = index_two(i, j, wf%n_ao)
           C_o(i, j) = wf%mo_coef(ij, 1) 
        enddo
!
        do j = 1, n_cc2_v 
           ij = index_two(i, j + wf%n_o, wf%n_ao)
           C_v(i, j) = wf%mo_coef(ij, 1) 
        enddo
!
     enddo
     write(unit_output,*)'10'
      flush(unit_output)
!
      call allocator(C_o_transformed, wf%n_ao, n_cc2_o)
      call dgemm('N', 'N',    &
                  wf%n_ao,    &
                  n_cc2_o,    &
                  n_cc2_o,    &
                  one,        &
                  C_o,        &
                  wf%n_ao,    &
                  M,          &
                  n_cc2_o,    &
                  zero,       &
                  C_o_transformed, &
                  wf%n_ao)
!
      call deallocator(C_o, wf%n_ao, n_cc2_o)
!
      call allocator(C_v_transformed, wf%n_ao, n_cc2_v)
      call dgemm('N', 'N',                   &
                  wf%n_ao,                   &
                  n_cc2_v,                   &
                  n_cc2_v,                   &
                  one,                       &
                  C_v,                       &
                  wf%n_ao,                   &
                  N,                         &
                  n_cc2_v,                   &
                  zero,                      &
                  C_v_transformed,           &
                  wf%n_ao)
!
      call deallocator(C_v, wf%n_ao, n_cc2_v)
      call deallocator(N, n_cc2_v, n_cc2_v)
      call deallocator(M, n_cc2_o, n_cc2_o)  
      write(unit_output,*)'11'
      flush(unit_output)
!
      do i = 1, wf%n_ao
!

         do j = 1, n_cc2_o
            ij = index_two(i, j, wf%n_ao)
            wf%mo_coef(ij, 1) = C_o_transformed(i, j) 
         enddo
!
         do j = 1, n_cc2_v 
            ij = index_two(i, j + wf%n_o, wf%n_ao)
            wf%mo_coef(ij, 1) = C_v_transformed(i, j) 
         enddo
!
      enddo
!
      write(unit_output,*)'12'
      flush(unit_output)
      call deallocator(C_o_transformed, wf%n_ao, n_cc2_o)
      call deallocator(C_v_transformed, wf%n_ao, n_cc2_v)
!
!     :: Determine number of active orbitals ::
!
      sum_o       = 1
      wf%n_CCSD_o = 1
!
      do while ((sum_o .gt. wf%mlcc_settings%delta_o) .and. (wf%n_CCSD_o .lt. n_cc2_o))
!
         sum_o = sum_o - eigenvalues_o(n_CC2_o - (wf%n_CCSD_o - 1), 1)
         wf%n_CCSD_o = wf%n_CCSD_o + 1
!
      enddo

!
      sum_v      = 1
      wf%n_CCSD_v = 1
!
      do while (sum_v .gt. wf%mlcc_settings%delta_v .and. (wf%n_CCSD_v .lt. n_cc2_v))
!
         sum_v = sum_v - eigenvalues_v(n_CC2_v - (wf%n_CCSD_v - 1), 1)
         wf%n_CCSD_v = wf%n_CCSD_v + 1
!
      enddo 
      write(unit_output,*)'13'
      flush(unit_output)
!
!     Save information to object
!
      wf%first_CCSD_o = 1
      wf%first_CCSD_v = 1
!
      wf%n_CC2_o = n_cc2_o - wf%n_CCSD_o
      wf%n_CC2_v = n_cc2_v - wf%n_CCSD_v
!
      wf%first_CC2_o = 1 + wf%n_CCSD_o
      wf%first_CC2_v = 1 + wf%n_CCSD_v
!
      call deallocator(eigenvalues_o, n_cc2_o, 1)
      call deallocator(eigenvalues_v, n_cc2_v, 1)
      call wf%print_cnto_info
!     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     -::- Finding orbital energies and new block diagonal C matrix -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     Initialize amplitudes and associated attributes
!
      call wf%initialize_amplitudes
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf%read_transform_cholesky
!
!     Initialize fock matrix
!
      call wf%initialize_fock_matrix
!
!     :: Occupied Orbitals ::
!
!     Diagonalize active-active block
!
      if (wf%n_CCSD_o .gt. 0) then
!
         call allocator(work, 4*(wf%n_CCSD_o), 1)
         call allocator(orbital_energies, (wf%n_CCSD_o), 1)
         work = zero
!
         call dsyev('V','U',              &
                     (wf%n_CCSD_o),       &
                     wf%fock_ij,          &
                     wf%n_o,              &
                     orbital_energies,    &
                     work,                & 
                     4*(wf%n_CCSD_o),     &
                     info)
!
         call deallocator(work, 4*(wf%n_CCSD_o), 1)
!
         if (info .ne. 0) then
            write(unit_output,*)'WARNING: Diagonalization of active virtual block not successful. '
            stop
         endif
!
!        Setting orbital energies
!
         do j = 1, wf%n_CCSD_o
!
            wf%fock_diagonal(j,1) = orbital_energies(j,1)
!
         enddo
!
         call deallocator(orbital_energies, (wf%n_CCSD_o), 1)
      endif
!
!     Diagonalize inactive-inactive block 
!
      if (wf%n_CC2_o .gt. 0) then
         call allocator(work, 4*wf%n_CC2_o, 1)
         call allocator(orbital_energies, wf%n_CC2_o, 1)
         orbital_energies = zero
         work = zero
!
         call dsyev('V','U',                                   &
                     wf%n_CC2_o,                               &
                     wf%fock_ij(wf%first_CC2_o, wf%first_CC2_o),&
                     wf%n_o,                                   &
                     orbital_energies,                         &
                     work,                                     & 
                     4*wf%n_CC2_o,                             &
                     info)
!
         call deallocator(work, 4*wf%n_CC2_o, 1)
!
         if (info .ne. 0) then
            write(unit_output,*)'WARNING: Diagonalization of inactive virtual block not successful.'
            stop
         endif
!
!        Setting orbital energies
!
         do j = 1, wf%n_CC2_o
!
            wf%fock_diagonal(j + wf%first_CC2_o - 1,1) = orbital_energies(j,1)
!
         enddo
         call deallocator(orbital_energies, wf%n_CC2_o, 1)
      endif
!
!     Transform C-matrix (active occupied block)
!
      call allocator(C_o, wf%n_ao, wf%n_o)
      C_o = zero     
!  
      do i = 1, wf%n_ao
!  
         do j = 1, wf%n_o
           ij = index_two(i, j, wf%n_ao)
           C_o(i, j) = wf%mo_coef(ij, 1) 
         enddo
!  
      enddo
!
      if (wf%n_CCSD_o .gt. 0) then
!  
         call allocator(C_o_transformed, wf%n_ao, wf%n_CCSD_o)
!  
         call dgemm('N', 'N',          &
                     wf%n_ao,          &
                     wf%n_CCSD_o,      &
                     wf%n_CCSD_o,      &
                     one,              &
                     C_o,              &
                     wf%n_ao,          &
                     wf%fock_ij,       &
                     wf%n_o,           &
                     zero,             &
                     C_o_transformed,  &
                     wf%n_ao)
!  
         do i = 1, wf%n_ao
!  
            do j = 1, wf%n_CCSD_o
               ij = index_two(i, j, wf%n_ao)
               wf%mo_coef(ij, 1) = C_o_transformed(i, j) 
            enddo
!  
         enddo
!  
         call deallocator(C_o_transformed, wf%n_ao, wf%n_CCSD_o)
      endif
!
!     Transform C-matrix (inactive occupied block)
!
      if (wf%n_CC2_o .gt. 0) then
         call allocator(C_o_transformed, wf%n_ao, wf%n_CC2_o)
         call dgemm('N', 'N',    &
                     wf%n_ao,    &
                     wf%n_CC2_o, &
                     wf%n_CC2_o, &
                     one,        &
                     C_o(1, wf%first_CC2_o),        &
                     wf%n_ao,    &
                     wf%fock_ij(wf%first_CC2_o, wf%first_CC2_o), &
                     wf%n_o,     &
                     zero,       &
                     C_o_transformed, &
                     wf%n_ao)
!
         call deallocator(C_o, wf%n_ao, wf%n_o)
!
         do i = 1, wf%n_ao
!
            do j = 1, wf%n_CC2_o
               ij = index_two(i, j + wf%first_CC2_o -1, wf%n_ao)
               wf%mo_coef(ij, 1) = C_o_transformed(i, j) 
            enddo
!
         enddo
!
         call deallocator(C_o_transformed, wf%n_ao, wf%n_CC2_o)
      endif
!
!     :: Vacant orbitals ::
!
!     Diagonalize active-active block
!
      if (wf%n_CCSD_v .gt. 0) then
         call allocator(work, 4*(wf%n_CCSD_v), 1)
         call allocator(orbital_energies, (wf%n_CCSD_v), 1)
         work = zero
!
         call dsyev('V','U',              &
                     (wf%n_CCSD_v),   &
                     wf%fock_ab,          &
                     wf%n_v,              &
                     orbital_energies,    &
                     work,                & 
                     4*(wf%n_CCSD_v), &
                     info)
!
         call deallocator(work, 4*(wf%n_CCSD_v), 1)
!
         if (info .ne. 0) then
            write(unit_output,*)'WARNING: Diagonalization of active virtual block not successful. '
                  stop
         endif
!
!        Setting orbital energies
!
         do j = 1, wf%n_CCSD_v
!
            wf%fock_diagonal(j + wf%n_o ,1) = orbital_energies(j,1)
!
         enddo
!
         call deallocator(orbital_energies, (wf%n_CCSD_v), 1)
!
      endif
!
!     Diagonalize inactive-inactive block 
!
      if (wf%n_CC2_v .gt. 0) then
         call allocator(work, 4*wf%n_CC2_v, 1)
         call allocator(orbital_energies, wf%n_CC2_v, 1)
         orbital_energies = zero
         work = zero
!
         call dsyev('V','U',                                      &
                     wf%n_CC2_v,                                  &
                     wf%fock_ab(wf%first_CC2_v, wf%first_CC2_v),  &
                     wf%n_v,                                      &
                     orbital_energies,                            &
                     work,                                        & 
                     4*(wf%n_CC2_v),                              &
                     info)
!
         call deallocator(work, 4*wf%n_CC2_v, 1)
!
         if (info .ne. 0) then
            write(unit_output,*)'WARNING: Diagonalization of inactive virtual block not successful.'
            stop
         endif
!
!        Setting orbital energies
!
         do j = 1, wf%n_CC2_v
!
            wf%fock_diagonal(j + wf%n_o + wf%first_CC2_v - 1,1) = orbital_energies(j,1)
!
         enddo

         call deallocator(orbital_energies, wf%n_CC2_v, 1)
!
      endif
!
!     Transform C-matrix (active virtual block)
!
      call allocator(C_v, wf%n_ao, wf%n_v)
      C_v = zero     
!
      do i = 1, wf%n_ao
!
         do j = 1, wf%n_v
            ij = index_two(i, j + wf%n_o, wf%n_ao)
            C_v(i, j) = wf%mo_coef(ij, 1) 
         enddo
!
      enddo
!
      if (wf%n_CCSD_v .gt. 0) then
!
         call allocator(C_v_transformed, wf%n_ao, wf%n_CCSD_v)
         call dgemm('N', 'N',    &
                     wf%n_ao,    &
                     wf%n_CCSD_v,&
                     wf%n_CCSD_v,&
                     one,        &
                     C_v,        &
                     wf%n_ao,    &
                     wf%fock_ab, &
                     wf%n_v,     &
                     zero,       &
                     C_v_transformed, &
                     wf%n_ao)
!
         do i = 1, wf%n_ao
            do j = 1, wf%n_CCSD_v
               ij = index_two(i, j + wf%n_o, wf%n_ao)
               wf%mo_coef(ij, 1) = C_v_transformed(i, j) 
            enddo
         enddo
!
         call deallocator(C_v_transformed, wf%n_ao, wf%n_CCSD_v)
!
      endif
!
!     Transform C-matrix (inactive virtual block)
!
      if (wf%n_CC2_v .gt. 0) then
!
         call allocator(C_v_transformed, wf%n_ao, wf%n_CC2_v)
         call dgemm('N', 'N',    &
                     wf%n_ao,    &
                     wf%n_CC2_v, &
                     wf%n_CC2_v, &
                     one,        &
                     C_v(1, wf%first_CC2_v),        &
                     wf%n_ao,    &
                     wf%fock_ab(wf%first_CC2_v, wf%first_CC2_v), &
                     wf%n_v,     &
                     zero,       &
                     C_v_transformed, &
                     wf%n_ao)
!
         do i = 1, wf%n_ao

            do j = 1, wf%n_CC2_v
               ij = index_two(i, j + wf%n_o + wf%first_CC2_v - 1, wf%n_ao)
               wf%mo_coef(ij, 1) = C_v_transformed(i, j) 
            enddo
         enddo
!
         call deallocator(C_v_transformed, wf%n_ao, wf%n_CC2_v)
!
      endif
!
      call deallocator(C_v, wf%n_ao, wf%n_v)
      
!
   end subroutine ccsd_cnto_cvs_mlccsd
!
!
   module subroutine print_cnto_info_mlccsd(wf)
!!
!!
      implicit none 
!
      class(mlccsd) :: wf
!
      write(unit_output,'(t3,a40, i3)') 'Number of CCSD occupied orbitals:      ', wf%n_CCSD_o
      write(unit_output,'(t3,a40, i3/)')'Number of CCSD virtual orbitals:       ', wf%n_CCSD_v
!
      write(unit_output,'(t3,a40, i3)') 'Number of CC2 occupied orbitals:       ', wf%n_CC2_o
      write(unit_output,'(t3,a40, i3/)')'Number of CC2 virtual orbitals:        ', wf%n_CC2_v
!
      write(unit_output,'(t3,a40, i3)') 'Number of inactive occupied orbitals:  ', wf%n_CCS_o
      write(unit_output,'(t3,a40, i3/)')'Number of inactive virtual orbitals:   ', wf%n_CCS_v
      flush(unit_output)
!
   end subroutine print_cnto_info_mlccsd
!
end submodule orbital_partitioning