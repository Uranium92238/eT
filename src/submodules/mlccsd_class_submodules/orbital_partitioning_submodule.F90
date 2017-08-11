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
   subroutine cholesky_localization_drv_mlccsd(wf)
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
   subroutine cholesky_localization_CCSD_CC2_CCS_mlccsd(wf, ao_center_info, n_ao_on_center,&
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
!     Looping variables
!
      integer(i15) :: i = 0
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
   subroutine cholesky_localization_CCSD_CCS_mlccsd(wf, ao_center_info, n_ao_on_center,&
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
       write(unit_output, *)'2'
      flush(unit_output)
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
   subroutine cholesky_localization_CCSD_CC2_mlccsd(wf, ao_center_info, n_ao_on_center,&
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
   subroutine construct_MO_transformation_matrix_mlccsd(wf)
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
end submodule orbital_partitioning