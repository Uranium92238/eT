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
   logical :: timings = .false.
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
!     Get number of active spaces
!
      wf%n_active_spaces = get_number_of_active_spaces(unit_cholesky_decomp)
!
!     Initialize orbital info variables
!
      call wf%initialize_orbital_info
!
!     Print number of active spaces requested active spaces
!
      write(unit_output,'(t3,a,i3,a/)')'Requested ', wf%n_active_spaces,' active space(s).'
!
!     Must keep information on CC2/CCS partition
!  
      call allocator(wf%mo_coef_cc2_ccs, wf%n_ao, wf%n_mo) 
      call allocator(wf%fock_diagonal_cc2_ccs, wf%n_mo, 1)
!
      call allocator(orbitals, wf%n_ao, wf%n_mo)
      call allocator(orbital_energies, wf%n_mo, 1)
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
         if (wf%n_active_spaces .gt. 1) then
            write(unit_output,*)'Not yet implemented'
            stop
         endif
!
         call wf%cholesky_localization_CCSD_CC2(ao_center_info, n_ao_on_center,&
                                                   ao_fock, n_nuclei, unit_cholesky_decomp)
!
      endif
!
!     CHECK ORTHOGONALITY!!
!
!
!     Deallocate ao-fock
!
      call deallocator(ao_fock, wf%n_ao, wf%n_ao)
      call deallocator_int(n_ao_on_center, n_nuclei, 1)      
      call deallocator_int(ao_center_info, wf%n_ao, 2)
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
      integer(i15) :: active_space = 0, i = 0
!
!     Active space variables
!
      integer(i15) :: n_active_aos
      integer(i15) :: n_active_orbitals_o
      integer(i15) :: n_active_orbitals_v
      integer(i15) :: n_CC2_atoms, n_CCSD_atoms
      integer(i15) :: n_vectors_o, n_vectors_v
      integer(i15) :: offset_o, offset_v
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
!     Start loop over active spaces
!  
      do active_space = 1, wf%n_active_spaces
!  
!        Get CC2/CCSD-active atoms
!  
         n_CC2_atoms  =  get_number_of_active_atoms(unit_cholesky_decomp, active_space, 'CC2  ')
         n_CCSD_atoms =  get_number_of_active_atoms(unit_cholesky_decomp, active_space, 'CCSD ')
!  
!        Allocate active atoms list
!  
         call allocator_int(active_atoms_CC2, n_CC2_atoms, 1)
         call allocator_int(active_atoms_CCSD, n_CCSD_atoms, 1)
         call allocator_int(active_atoms, n_CCSD_atoms + n_CC2_atoms, 1)
!  
!        Get list of active atoms
! 
         call get_active_atoms(unit_cholesky_decomp, active_atoms_CC2,  n_CC2_atoms,  active_space, 'CC2  ')
         call get_active_atoms(unit_cholesky_decomp, active_atoms_CCSD, n_CCSD_atoms, active_space, 'CCSD ')
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
!        Sanity check on active atoms
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
!        :: CC2 localized Cholesky orbitals ::
!  
         n_active_aos = 0
!
         do i = 1, n_CC2_atoms + n_CCSD_atoms
            n_active_aos = n_active_aos + n_ao_on_center(active_atoms(i,1),1)
         enddo
!  
!        Construct active_ao_index_list
!  
         call allocator_int(active_ao_index_list, n_active_aos, 1)
         call construct_active_ao_index_list(active_ao_index_list, n_active_aos, active_atoms, &
                                                  n_CC2_atoms + n_CCSD_atoms, ao_center_info, wf%n_ao)
!  
!        Occupied part
!
         n_vectors_o = 0  
         call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                              density_o, n_vectors_o,&
                                             .true., n_active_aos, active_ao_index_list)
!  
!        Virtual part
!  
         n_vectors_v = 0 
         call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                              density_v, n_vectors_v,&
                                              .true., n_active_aos, active_ao_index_list)
!  
!        Calculate new offset         
!  
         offset_o = offset_o + n_vectors_o
         offset_v = offset_v + n_vectors_v
!  
         call deallocator_int(active_atoms, n_CC2_atoms + n_CCSD_atoms, 1)
         call deallocator_int(active_ao_index_list, n_active_aos, 1)
!     
      enddo
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
      do active_space = 1, wf%n_active_spaces
!
!        Construct CC2 density matrix, occupied and virtual
!
         call allocator(density_v, wf%n_ao, wf%n_ao)   
         call allocator(density_o, wf%n_ao, wf%n_ao)
!
!        Construct CC2 density
!
         call allocator(C_cc2, wf%n_ao, wf%n_CC2_o(active_space,1) + wf%n_CC2_v(active_space,1)) 
!
         do i = 1, wf%n_CC2_o(active_space,1)
            C_cc2(:,i) = wf%mo_coef_cc2_ccs(:,i)
         enddo
!
         do i = 1, wf%n_CC2_v(active_space,1)
            C_cc2(:, i + wf%n_CC2_o(active_space,1)) = wf%mo_coef_cc2_ccs(:, i + wf%n_o)
         enddo
!
         call wf%construct_density_matrices(density_o, density_v, C_cc2, &
                                       wf%n_CC2_o(active_space,1), wf%n_CC2_v(active_space,1))
!  
         call deallocator(C_cc2, wf%n_ao, wf%n_CC2_o(active_space,1) + wf%n_CC2_v(active_space,1)) 
!
         n_CCSD_atoms =  get_number_of_active_atoms(unit_cholesky_decomp, active_space, 'CCSD ')
!  
!        Allocate active atoms list
!
         call allocator_int(active_atoms, n_CCSD_atoms, 1)
!  
!        Get list of active atoms
! 
         call get_active_atoms(unit_cholesky_decomp, active_atoms,  n_CCSD_atoms, active_space, 'CCSD ')
!  
!        Construct active_ao_index_list
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
!        Occupied part
!  
         n_vectors_o = 0 
         call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                              density_o, n_vectors_o,&
                                             .true., n_active_aos, active_ao_index_list)
!  
!        Virtual part
!
         n_vectors_v = 0 
         call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                              density_v, n_vectors_v,&
                                              .true., n_active_aos, active_ao_index_list)
!
         call deallocator_int(active_ao_index_list, n_active_aos, 1)
!
!  
!        Save active space information
!  
         wf%n_CCSD_o(active_space, 1) = n_vectors_o
         wf%n_CCSD_v(active_space, 1) = n_vectors_v
!
         wf%first_CCSD_o(active_space, 1) = offset_o
         wf%first_CCSD_v(active_space, 1) = offset_v - wf%n_o
!  
!        Calculate new offset         
!  
         offset_o = offset_o + n_vectors_o
         offset_v = offset_v + n_vectors_v  
!
!        :: CC2  localized Cholesky orbitals  ::
!  
!        Occupied part
!  
         n_vectors_o = 0 
         call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                              density_o, n_vectors_o,&
                                             .false., n_active_aos)
!  
!        Virtual part
!  
         n_vectors_v = 0 
         call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                              density_v, n_vectors_v,&
                                              .false., n_active_aos)
!
!        Save inactive space information
!
         wf%n_CC2_o(active_space,1) = n_vectors_o
         wf%n_CC2_v(active_space,1) = n_vectors_v
!
         wf%first_CC2_o(active_space,1) = offset_o
         wf%first_CC2_v(active_space,1) = offset_v - wf%n_o
!  
!        Calculate new offset         
!  
         offset_o = offset_o + n_vectors_o
         offset_v = offset_v + n_vectors_v
!
         call deallocator(density_v, wf%n_ao, wf%n_ao)   
         call deallocator(density_o, wf%n_ao, wf%n_ao)
!
      enddo
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
      integer(i15) :: active_space = 0, i = 0
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
      call allocator (orbitals, wf%n_ao, wf%n_mo)
      call allocator (orbital_energies, wf%n_mo, 1)
!
      wf%mo_coef_cc2_ccs       = wf%mo_coef
      wf%fock_diagonal_cc2_ccs = wf%fock_diagonal
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
!     Start loop over active spaces
!  
      do active_space = 1, wf%n_active_spaces
!  
!        Get CC2/CCSD-active atoms
!
         n_CCSD_atoms =  get_number_of_active_atoms(unit_cholesky_decomp, active_space, 'CCSD ')  
!
!        Allocate active atoms list
!  
         call allocator_int(active_atoms, n_CCSD_atoms, 1)
!  
!        Get list of active atoms
!  
         call get_active_atoms(unit_cholesky_decomp, active_atoms,  n_CCSD_atoms, active_space, 'CCSD  ')
!  
!        Sanity check on active atoms
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
!        :: Constructing active (CCSD) localized Cholesky orbitals ::
!  
         n_active_aos = 0
         do i = 1, n_CCSD_atoms
            n_active_aos = n_active_aos + n_ao_on_center(active_atoms(i,1),1)
         enddo
!  
!        Construct active_ao_index_list
!  
         call allocator_int(active_ao_index_list, n_active_aos, 1)
         call construct_active_ao_index_list(active_ao_index_list, n_active_aos, active_atoms, &
                                                  n_CCSD_atoms, ao_center_info, wf%n_ao)
!  
!        Occupied part
!  
         n_vectors_o = 0 
         call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                              density_o, n_vectors_o,&
                                             .true., n_active_aos, active_ao_index_list)
!  
!        Virtual part
!  
         n_vectors_v = 0 
         call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                              density_v, n_vectors_v,&
                                              .true., n_active_aos, active_ao_index_list)
!  
!        Save active space information
!  
         wf%n_CCSD_o(active_space, 1) = n_vectors_o
         wf%n_CCSD_v(active_space, 1) = n_vectors_v
!
         wf%first_CCSD_o(active_space, 1) = offset_o
         wf%first_CCSD_v(active_space, 1) = offset_v - wf%n_o
!  
!        Calculate new offset         
!  
         offset_o = offset_o + n_vectors_o
         offset_v = offset_v + n_vectors_v
!  
         call deallocator_int(active_atoms, n_CCSD_atoms, 1)
         call deallocator_int(active_ao_index_list, n_active_aos, 1)
!     
      enddo
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
      integer(i15) :: active_space = 0, i = 0
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
        call allocator (orbitals, wf%n_ao, wf%n_mo)
         call allocator (orbital_energies, wf%n_mo, 1)
!
         wf%mo_coef_cc2_ccs       = wf%mo_coef
         wf%fock_diagonal_cc2_ccs = wf%fock_diagonal
!
!        ::::::::::::::::::::::::::::::::
!        -::-   CCSD/CC2 partition   -::-
!        ::::::::::::::::::::::::::::::::
!
!        :: Construct canonical occupied and vacant density matrices ::    
!
         call allocator(density_o, wf%n_ao, wf%n_ao)
         call allocator(density_v, wf%n_ao, wf%n_ao)
         density_o = zero
         density_v = zero
!
         call wf%construct_density_matrices(density_o, density_v, wf%mo_coef, wf%n_o, wf%n_v)
!
!        No CCS region
!
         wf%n_CCS_o = 0
         wf%n_CCS_v = 0
!
         offset_o = 1
         offset_v = 1 + wf%n_o
!
!        Start loop over active spaces
!  
         do active_space = 1, wf%n_active_spaces
!  
!           Get CC2/CCSD-active atoms
!
            n_CCSD_atoms =  get_number_of_active_atoms(unit_cholesky_decomp, active_space, 'CCSD ')  
!
!           Allocate active atoms list
!  
            call allocator_int(active_atoms, n_CCSD_atoms, 1)
!  
!           Get list of active atoms
!  
            call get_active_atoms(unit_cholesky_decomp, active_atoms,  n_CCSD_atoms, active_space, 'CCSD  ')
!  
!           Sanity check on active atoms
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
!           :: Constructing active (CCSD) localized Cholesky orbitals ::
!  
            n_active_aos = 0
            do i = 1, n_CCSD_atoms
               n_active_aos = n_active_aos + n_ao_on_center(active_atoms(i,1),1)
            enddo
!  
!           Construct active_ao_index_list
!  
            call allocator_int(active_ao_index_list, n_active_aos, 1)
            call construct_active_ao_index_list(active_ao_index_list, n_active_aos, active_atoms, &
                                                     n_CCSD_atoms, ao_center_info, wf%n_ao)
!  
!           Occupied part
!  
            n_vectors_o = 0 
            call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                                 density_o, n_vectors_o,&
                                                .true., n_active_aos, active_ao_index_list)
!  
!           Virtual part
!  
            n_vectors_v = 0 
            call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                                 density_v, n_vectors_v,&
                                                 .true., n_active_aos, active_ao_index_list)
!  
!           Save active space information
!  
            wf%n_CCSD_o(active_space, 1) = n_vectors_o
            wf%n_CCSD_v(active_space, 1) = n_vectors_v
!
            wf%first_CCSD_o(active_space, 1) = offset_o
            wf%first_CCSD_v(active_space, 1) = offset_v - wf%n_v
!  
!           Calculate new offset         
!  
            offset_o = offset_o + n_vectors_o
            offset_v = offset_v + n_vectors_v
!  
            call deallocator_int(active_atoms, n_CCSD_atoms, 1)
            call deallocator_int(active_ao_index_list, n_active_aos, 1)
!     
         enddo
!
!        :: CC2  localized Cholesky orbitals  ::
!
!        Occupied part
!
         n_vectors_o = 0 
         call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                           density_o, n_vectors_o,&
                                           .false., n_active_aos)
!
!        Virtual part
!
         n_vectors_v = 0 
         call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                           density_v, n_vectors_v,&
                                           .false., n_active_aos)
!
!        Save inactive space information
!
         wf%n_CC2_o(1,1) = n_vectors_o
         wf%n_CC2_v(1,1) = n_vectors_v
!
         wf%first_CC2_o(1,1) = offset_o
         wf%first_CC2_v(1,1) = offset_v - wf%n_o
!
         wf%mo_coef       = orbitals
         wf%fock_diagonal = orbital_energies
!
!        Deallocations
!
         call deallocator(density_v, wf%n_ao, wf%n_ao)   
         call deallocator(density_o, wf%n_ao, wf%n_ao)
!
   end subroutine cholesky_localization_CCSD_CC2_mlccsd
!
!
   subroutine construct_MO_transformation_matrix(wf)
!!
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(:,:), allocatable :: X
      real(dp), dimension(:,:), allocatable :: S_packed, S
!
      integer(i15) :: unit_overlap = -1, ioerror = 0
!
      call allocator(wf%T, wf%n_total_active_v + wf%n_total_active_o, wf%n_total_active_o + wf%n_total_active_v)
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
      call allocator(X, wf%n_total_active_v+wf%n_total_active_o, wf%n_ao)
!
      call dgemm('T', 'N',                                 &
                  wf%n_total_active_v+wf%n_total_active_o, &
                  wf%n_ao,                                 &
                  wf%n_ao,                                 &
                  one,                                     &
                  wf%mo_coef,                              &
                  wf%n_ao,                                 &
                  S,                                       &
                  wf%n_ao,                                 &
                  zero,                                    &
                  X,                                       &
                  wf%n_total_active_v+wf%n_total_active_o)
!
      call deallocator(S, wf%n_ao, wf%n_ao)
!
      call dgemm('N', 'N',                                 &
                  wf%n_total_active_v+wf%n_total_active_o, &
                  wf%n_total_active_v+wf%n_total_active_o, &
                  wf%n_ao,                                 &
                  one,                                     &
                  X,                                       &
                  wf%n_mo,                                 &
                  wf%mo_coef_cc2_ccs,                      &
                  wf%n_total_active_v+wf%n_total_active_o, &
                  zero,                                    &
                  wf%T,                                    &
                  wf%n_total_active_v+wf%n_total_active_o)
!
      call deallocator(X, wf%n_mo, wf%n_ao)
!
      call deallocator(wf%mo_coef_cc2_ccs, wf%n_ao, wf%n_mo)
!
   end subroutine construct_MO_transformation_matrix
!
!
end submodule orbital_partitioning