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
   module subroutine orbital_partitioning_mlccsd(wf)
!!
!!    Orbital partitioning,
!!    Written by Sarai D. Folkestad, June 2017
!!
!!    Directs the partitioning for mlcc calculations.
!!
      implicit none
!
      class(mlccsd) :: wf
!
      if (wf%CCSD_orbitals%cholesky) then
!
!        If cholesky - do Cholesky decomposition
!
         call wf%cholesky_localization_drv
!
      elseif (wf%CCSD_orbitals%cnto) then
!
!        If CNTO - do CNTO
!     
         if (.not. wf%mlcc_settings%CC2) then
!
            write(unit_output,*)'Error: CC2 space required for MLCCSD with cntos'
            stop
!
         endif
!
         call wf%cnto_orbital_drv
!
      endif
!
   end subroutine orbital_partitioning_mlccsd
!
!
   module subroutine cholesky_localization_drv_mlccsd(wf)
!!
!!    Cholesky orbital localization driver,
!!    Written by Sarai D. Folkestad, July 2017.
!!
!!    Driver for Cholesky density decomposition 
!!
!!    - Collects atom and ao-basis information.
!!    - Constructs occupied and vacant densities.
!!    - Constructs AO Fock matrix.  (This is currently an N^5 operation, should be optimized/removed)
!!    - By looping over active spaces, the occupied and virtual densities are Cholesky decomposed
!!      and the cholesky vectors are used to generate new localized MO's.
!!    - New orbitals are tested for orthonormality (Not implemented yet, only need overlap matrix from DALTON)    
!!
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
      write(unit_output,'(//t3,a)')    ':: Cholesky orbital construction'
      write(unit_output,'(t3,a/)')   ':: E. F. Kjønstad, S. D. Folkestad, Jun 2017'
      flush(unit_output)
!
!     Get center info
!
      call read_atom_info(n_nuclei, wf%n_ao)
!
      call wf%mem%alloc_int(n_ao_on_center, n_nuclei, 1)      
      call wf%mem%alloc_int(ao_center_info, wf%n_ao, 2)
!
      call read_center_info(n_nuclei, wf%n_ao, n_ao_on_center, ao_center_info)
!
!     Start timings
!
      call cpu_time(start_chol_deco)
!
!     Construct AO fock matrix
!
      call wf%mem%alloc(density_o, wf%n_ao, wf%n_ao)
      call wf%mem%alloc(density_v, wf%n_ao, wf%n_ao)
      density_o = zero
      density_v = zero
!
      call wf%construct_density_matrices(density_o, density_v, wf%mo_coef, wf%n_o, wf%n_v)
!
      call wf%mem%alloc(ao_fock, wf%n_ao, wf%n_ao)
      call wf%construct_ao_fock(ao_fock)
!
      call wf%mem%dealloc(density_o, wf%n_ao, wf%n_ao)
      call wf%mem%dealloc(density_v, wf%n_ao, wf%n_ao)
!
!     Test for CCS region
!    
      if (wf%mlcc_settings%CCS) then
!
!        Test for CC2 region
!
         if (wf%mlcc_settings%CC2) then
                        call wf%cholesky_localization_CCSD_CC2_CCS(ao_center_info, n_ao_on_center,&
                                                       ao_fock, n_nuclei, unit_cholesky_decomp)
!
         else ! No CC2 region
!
            call wf%cholesky_localization_CCSD_CCS(ao_center_info, n_ao_on_center,&
                                                       ao_fock, n_nuclei, unit_cholesky_decomp)
!
         endif
!
      else ! No CCS region 
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
      call wf%mem%alloc(S_packed, wf%n_ao*(wf%n_ao+1)/2, 1)
      read(unit_overlap,*)S_packed
      close(unit_overlap)
!
      call wf%mem%alloc(S, wf%n_ao, wf%n_ao)
!
      call squareup(S_packed, S, wf%n_ao)
!
      call wf%mem%dealloc(S_packed, wf%n_ao*(wf%n_ao+1)/2, 1)
!
      call wf%mem%alloc(Y_1, wf%n_mo, wf%n_ao)
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
      call wf%mem%dealloc(S, wf%n_ao, wf%n_ao)
      call wf%mem%alloc(Y_2, wf%n_mo, wf%n_mo)
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
      call wf%mem%dealloc(Y_1, wf%n_mo, wf%n_ao)
!
      if (.not. (check_orthogonality(Y_2, wf%n_mo, wf%n_mo))) then
         write(unit_output,*)'New orbitals not orthogonal'
         stop
      endif
!
      call wf%mem%dealloc(Y_2, wf%n_mo, wf%n_mo)
!
!     Deallocate ao-fock
!
      call wf%mem%dealloc(ao_fock, wf%n_ao, wf%n_ao)
      call wf%mem%dealloc_int(n_ao_on_center, n_nuclei, 1)      
      call wf%mem%dealloc_int(ao_center_info, wf%n_ao, 2)
!
!     Print decomposition info
!
      write(unit_output,'(t3,a,a,a/)')'Summary of ', trim(wf%name), ' Cholesky orbital construction:'
!
!     Print timings
!
      call cpu_time(end_chol_deco)
      write(unit_output,'(t6,a25,f14.8/)') 'Total CPU time (seconds):',&
                            end_chol_deco - start_chol_deco
!
      call wf%print_orbital_info
!
      flush(unit_output)
!
   end subroutine cholesky_localization_drv_mlccsd
!
   module subroutine cholesky_localization_CCSD_CC2_CCS_mlccsd(wf, ao_center_info, n_ao_on_center,&
                                                       ao_fock, n_nuclei, unit_cholesky_decomp)
!!
!!    Cholesky orbital localization CCS/CC2/CCSD,
!!    Written by Sarai D. Folkestad, July 2017
!!
!!    Cholesky partitiining routine for CCS/CC2/CCSD calculation
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
      integer(i15), dimension(:,:), allocatable :: active_atoms
      integer(i15), dimension(:,:), allocatable :: active_ao_index_list
!
!     Active space variables
!
      integer(i15) :: n_active_aos
      integer(i15) :: n_active_orbitals_o
      integer(i15) :: n_active_orbitals_v
      integer(i15) :: n_vectors_o, n_vectors_v
      integer(i15) :: offset_o, offset_v, offset
!
!     :::::::::::::::::::::::::::::::::::::
!     -::-  Prepare for decomposition  -::-
!     :::::::::::::::::::::::::::::::::::::
!
      call wf%mem%alloc(orbitals, wf%n_ao, wf%n_mo)
      call wf%mem%alloc(orbital_energies, wf%n_mo, 1)
!
!     ::::::::::::::::::::::::::::::::
!     -::-    CC2/CCS partition   -::-
!     ::::::::::::::::::::::::::::::::
!
!     :: Construct canonical occupied and vacant density matrices ::    
!
      call wf%mem%alloc(density_o, wf%n_ao, wf%n_ao)
      call wf%mem%alloc(density_v, wf%n_ao, wf%n_ao)
      density_o = zero
      density_v = zero
!
      call wf%construct_density_matrices(density_o, density_v, wf%mo_coef, wf%n_o, wf%n_v)
!
      offset_o = 1
      offset_v = 1 + wf%n_o
!  
      call wf%mem%alloc_int(active_atoms, &
                      wf%CCSD_orbitals%n_active_atoms + wf%CC2_orbitals%n_active_atoms, 1)
!
      do i = 1, wf%CC2_orbitals%n_active_atoms
         active_atoms(i, 1) = wf%CC2_orbitals%active_atoms(i, 1)
      enddo
!
      do i = 1, wf%CCSD_orbitals%n_active_atoms
         active_atoms(wf%CC2_orbitals%n_active_atoms + i, 1) = wf%CCSD_orbitals%active_atoms(i,1)
      enddo
! 
!     Sanity check on active atoms
! 
       if ((wf%CC2_orbitals%n_active_atoms + wf%CCSD_orbitals%n_active_atoms) .gt. n_nuclei) then
         write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
         stop
      endif
! 
      do i = 1, wf%CC2_orbitals%n_active_atoms + wf%CCSD_orbitals%n_active_atoms
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
      do i = 1, wf%CC2_orbitals%n_active_atoms + wf%CCSD_orbitals%n_active_atoms
         n_active_aos = n_active_aos + n_ao_on_center(active_atoms(i,1),1)
      enddo
!  
!     Construct active_ao_index_list
!  
      call wf%mem%alloc_int(active_ao_index_list, n_active_aos, 1)
      call construct_active_ao_index_list(active_ao_index_list, n_active_aos, active_atoms, &
                                               wf%CC2_orbitals%n_active_atoms + wf%CCSD_orbitals%n_active_atoms,&
                                                ao_center_info, wf%n_ao)
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
      call wf%mem%dealloc_int(active_atoms, wf%CC2_orbitals%n_active_atoms + wf%CCSD_orbitals%n_active_atoms, 1)
      call wf%mem%dealloc_int(active_ao_index_list, n_active_aos, 1)   
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
!     Total number of active orbitals
!
      n_active_orbitals_o = wf%n_o - wf%n_CCS_o
      n_active_orbitals_v = wf%n_v - wf%n_CCS_v
!
      wf%mo_coef_cc2_ccs       = orbitals
      wf%fock_diagonal_cc2_ccs = orbital_energies
!
!     Deallocations
!
      call wf%mem%dealloc(density_v, wf%n_ao, wf%n_ao)   
      call wf%mem%dealloc(density_o, wf%n_ao, wf%n_ao)
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
      call wf%mem%alloc(density_v, wf%n_ao, wf%n_ao)   
      call wf%mem%alloc(density_o, wf%n_ao, wf%n_ao)
!
!     Construct CC2 density
!
      call wf%mem%alloc(C_cc2, wf%n_ao*(wf%n_CC2_o + wf%n_CC2_v), 1) 
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
      call wf%mem%dealloc(C_cc2, wf%n_ao*(wf%n_CC2_o + wf%n_CC2_v), 1) 
! 
!     Construct active_ao_index_list
!  
      n_active_aos = 0
!
      do i = 1, wf%CCSD_orbitals%n_active_atoms
         n_active_aos = n_active_aos + n_ao_on_center(wf%CCSD_orbitals%active_atoms(i,1),1)
      enddo
!
      call wf%mem%alloc_int(active_ao_index_list, n_active_aos, 1)
      active_ao_index_list = 0
!
      call construct_active_ao_index_list(active_ao_index_list, n_active_aos, wf%CCSD_orbitals%active_atoms, &
                                               wf%CCSD_orbitals%n_active_atoms, ao_center_info, wf%n_ao)
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
      call wf%mem%dealloc_int(active_ao_index_list, n_active_aos, 1)
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
      call wf%mem%dealloc(density_v, wf%n_ao, wf%n_ao)   
      call wf%mem%dealloc(density_o, wf%n_ao, wf%n_ao)
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
!!    Cholesky orbital localization CCS/CCSD,
!!    Written by Sarai D. Folkestad, July 2017
!!
!!    Cholesky partitiining routine for CCS/CCSD calculation
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
      integer(i15) :: n_vectors_o, n_vectors_v
      integer(i15) :: offset_o, offset_v
!
!     :::::::::::::::::::::::::::::
!     -::-  CCS/CCSD partition -::-
!     :::::::::::::::::::::::::::::
!
!
      call wf%mem%alloc(orbitals, wf%n_ao, wf%n_mo)
      call wf%mem%alloc(orbital_energies, wf%n_mo, 1)
!
!     :: Construct canonical occupied and vacant density matrices ::    
!
      call wf%mem%alloc(density_o, wf%n_ao, wf%n_ao)
      call wf%mem%alloc(density_v, wf%n_ao, wf%n_ao)
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
!     Sanity check on active atoms
!  
      if (wf%CCSD_orbitals%n_active_atoms .gt. n_nuclei) then
         write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
         stop
      endif
!  
      do i = 1, wf%CCSD_orbitals%n_active_atoms
         if (wf%CCSD_orbitals%active_atoms(i,1) .gt. n_nuclei) then
            write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
            stop
         endif
      enddo
!  
!     :: Constructing active (CCSD) localized Cholesky orbitals ::
!  
      n_active_aos = 0
      do i = 1, wf%CCSD_orbitals%n_active_atoms
         n_active_aos = n_active_aos + n_ao_on_center(wf%CCSD_orbitals%active_atoms(i,1),1)
      enddo
!  
!     Construct active_ao_index_list
!  
      call wf%mem%alloc_int(active_ao_index_list, n_active_aos, 1)
      call construct_active_ao_index_list(active_ao_index_list, n_active_aos, wf%CCSD_orbitals%active_atoms, &
                                               wf%CCSD_orbitals%n_active_atoms, ao_center_info, wf%n_ao)
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
      call wf%mem%dealloc_int(active_ao_index_list, n_active_aos, 1)
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
      call wf%mem%dealloc(density_v, wf%n_ao, wf%n_ao)   
      call wf%mem%dealloc(density_o, wf%n_ao, wf%n_ao)
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
!!    Cholesky orbital localization CC2/CCSD
!!    Written by Sarai D. Folkestad, July 2017
!!
!!    Cholesky partitiining routine for CC2/CCSD calculation
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
      integer(i15), dimension(:,:), allocatable :: active_ao_index_list
!
!     Looping variables
!
      integer(i15) :: i = 0, j = 0, ij = 0
!
      integer(i15) :: n_active_aos
      integer(i15) :: n_active_orbitals_o
      integer(i15) :: n_active_orbitals_v
      integer(i15) :: n_vectors_o, n_vectors_v
      integer(i15) :: offset_o, offset_v
!
      call wf%mem%alloc(orbitals, wf%n_ao, wf%n_mo)
      call wf%mem%alloc(orbital_energies, wf%n_mo, 1)
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
      call wf%mem%alloc(density_o, wf%n_ao, wf%n_ao)
      call wf%mem%alloc(density_v, wf%n_ao, wf%n_ao)
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
!     Sanity check on active atoms
!  
      if (wf%CCSD_orbitals%n_active_atoms .gt. n_nuclei) then
         write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
         stop
      endif
!  
      do i = 1, wf%CCSD_orbitals%n_active_atoms
         if (wf%CCSD_orbitals%active_atoms(i,1) .gt. n_nuclei) then
            write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
            stop
         endif
      enddo 
!  
!     :: Constructing active (CCSD) localized Cholesky orbitals ::
!  
      n_active_aos = 0
      do i = 1, wf%CCSD_orbitals%n_active_atoms
         n_active_aos = n_active_aos + n_ao_on_center(wf%CCSD_orbitals%active_atoms(i,1),1)
      enddo
!  
!     Construct active_ao_index_list
!  
      call wf%mem%alloc_int(active_ao_index_list, n_active_aos, 1)
      call construct_active_ao_index_list(active_ao_index_list, n_active_aos, wf%CCSD_orbitals%active_atoms, &
                                                  wf%CCSD_orbitals%n_active_atoms, ao_center_info, wf%n_ao) 
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
      call wf%mem%dealloc_int(active_ao_index_list, n_active_aos, 1)
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
      call wf%mem%dealloc(density_v, wf%n_ao, wf%n_ao)   
      call wf%mem%dealloc(density_o, wf%n_ao, wf%n_ao)
!
      call wf%mem%dealloc(orbitals, wf%n_ao, wf%n_mo)
      call wf%mem%dealloc(orbital_energies, wf%n_mo, 1) 
!
   end subroutine cholesky_localization_CCSD_CC2_mlccsd
!
!
   module subroutine construct_MO_transformation_matrix_mlccsd(wf)
!!
!!       Construct MO transformation matrix,
!!       Written by Sarai D. Fokestad, July 2017
!!
!!       Constructs transformation matrix,
!!
!!          T = (C_CCSD)^T * S * C_CC2
!!
!!       between CC2 basis and CCSD basis. 
!!       Needed for transforming s_ij_ab from CC2 to CCSD basis.
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
      call wf%mem%alloc(S_packed, wf%n_ao*(wf%n_ao+1)/2,1)
      read(unit_overlap,*)S_packed
      close(unit_overlap)
!
      call wf%mem%alloc(S, wf%n_ao, wf%n_ao)
!
      call squareup(S_packed, S, wf%n_ao)
!
      call wf%mem%dealloc(S_packed, wf%n_ao*(wf%n_ao+1)/2,1)
!
      call wf%mem%alloc(X, wf%n_mo, wf%n_ao)
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
      call wf%mem%dealloc(S, wf%n_ao, wf%n_ao)
      call wf%mem%alloc(X2, wf%n_mo, wf%n_mo)
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
      call wf%mem%dealloc(X, wf%n_mo, wf%n_ao)
!
      call wf%mem%dealloc(wf%mo_coef_cc2_ccs, wf%n_ao, wf%n_mo)
      call wf%mem%alloc(wf%T_o, wf%n_o, wf%n_o)
      call wf%mem%alloc(wf%T_v, wf%n_v, wf%n_v)
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
      call wf%mem%dealloc(X2, wf%n_mo, wf%n_mo)
!
   end subroutine construct_MO_transformation_matrix_mlccsd
!
!
   module subroutine cnto_orbital_drv_mlccsd(wf)
!!
!!    CNTO orbital driver,
!!    Written by Sarai D. Folkestad, June 2017.
!!
!!    A CC2 calculation ground state and excited states is performed.
!!    The M and N matrices are then constructed, 
!! 
!!       M_i_j = sum_{a} R_ai R_aj + 1/2sum_{abk}(1+ δ_ai,bk δ_i,j )R_aibkR_ajbk
!!       N_a_b = sum_{i} R_ai R_bi + 1/2sum_{ijc}(1+ δ_ai,cj δ_a,b )R_aicjR_bicj.
!!   
!!    where R_ai is the i'th single part of the excitation vector obtained from the CC2 calculation
!!    and R_aibj is the doubles part of the excitation vectors.
!!
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
!!    the orbitals and orbital energies used in the CCSD calculation.
!!
      implicit none 
!
      class(mlccsd) :: wf
!
      real(dp) :: start_cnto = 0, end_cnto = 0
!
      integer(i15) :: cc2_n_parameters = 0, cc2_n_x2am = 0, n_cc2_o = 0, n_cc2_v = 0
!
!     Timings 
!
      call cpu_time(start_cnto)
!     
!     Prints
!     
      write(unit_output,'(//t3,a)') ':: CNTO orbital partitioning for MLCCSD calculation '
      write(unit_output,'(t3,a/)')   ':: E. F. Kjønstad, S. D. Folkestad, Jun 2017'
      write(unit_output,'(/a35,e10.2)')'Threshold for occupied orbitals:',wf%CCSD_orbitals%delta_o
      write(unit_output,'(a35,e10.2/)')'Threshold for virtual orbitals: ',wf%CCSD_orbitals%delta_v
!
!     :: Run lower level method ::
!
      if (wf%tasks%excited_state .or. wf%tasks%core_excited_state) then
!
         call wf%cnto_lower_level_method
!
      else
!
         write(unit_output,*)'Error: cntos without excited state calculation makes no sense.'
         stop
!
      endif
!
!     :: Construct orbitals ::
!
      call wf%cnto_orbitals
!
!     Print summary
!     
      write(unit_output,'(/t3,a,a,a/)')'Summary of ', trim(wf%name), ' CNTO orbital construction:'
!
!     Print timings
!
      call cpu_time(end_cnto)
!
      if (timings) write(unit_output,'(/t6,a26,f14.8/)') 'Total CPU time (seconds): ', end_cnto - start_cnto
!
      call wf%print_orbital_info 
      flush(unit_output)
!
   end subroutine cnto_orbital_drv_mlccsd
!
!
      module subroutine cnto_lower_level_method_mlccsd(wf)
!!
!!    CNTO lower level calculation (MLCCSD),
!!    Written by Sarai D. Folkestad, Aug. 2017
!!
!!    Runs lower level method for CNTOs
!!
      implicit none 
!
      class(mlccsd) :: wf
!
      type(mlcc2), allocatable :: cc2_wf
!
      integer(i15) :: i = 0, j = 0, ij = 0
!
!     ::::::::::::::::::::::::::::::::::::::::::::::::
!     -::- Running lower level method calculation -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::::
!
!     Allocate lower level method
!
      allocate(cc2_wf)
!
!     Initialize cc2_wf
!
      call wf%cnto_init_cc2(cc2_wf)
!
!     Call driver of lower level method
!
      call cc2_wf%drv
!
!     Set cc2 mo coefficients
!
      do i = 1, wf%n_ao
         do j = 1, wf%n_mo
            ij = index_two(i, j, wf%n_ao)
            wf%mo_coef_cc2_ccs(i, j) = cc2_wf%mo_coef(ij, 1)
         enddo
      enddo
!
!     Set cc2 fock matrix diagonal
!
      wf%fock_diagonal_cc2_ccs = cc2_wf%fock_diagonal
!
      wf%mo_coef = cc2_wf%mo_coef
      wf%fock_diagonal = cc2_wf%fock_diagonal
!
      wf%n_CC2_o = cc2_wf%n_CC2_o
      wf%n_CC2_v = cc2_wf%n_CC2_v
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
   end subroutine cnto_lower_level_method_mlccsd
!
!
   module subroutine cnto_orbitals_mlccsd(wf)
!!
!!    CNTO Oritals (MLCCSD),
!!    Written by Sarai D. Folkestad Aug. 2017
!!
!!    Constructs the CNTO orbitals based on exitation vectors from lower level method
!!
      implicit none
!
      class(mlccsd) :: wf
!
      integer(i15) :: cc2_n_x2am, cc2_n_parameters
      integer(i15) :: n_CC2_o, n_CC2_v
!
      integer(i15) :: unit_solution = -1
      integer(i15) :: ioerror = 0
      integer(i15) :: state
!
      real(dp), dimension(:,:), allocatable :: solution
      real(dp), dimension(:,:), allocatable :: R
      real(dp), dimension(:,:), allocatable :: R_restricted
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
!
      integer(i15) :: info
      integer(i15) :: i = 0, j = 0, a = 0, b = 0, c = 0, k = 0
      integer(i15) :: ai = 0, ij = 0, cj = 0, bj = 0, bk = 0, ai_full
      integer(i15) :: aibj = 0, aibk = 0, aib = 0, icj = 0
!
      real(dp) :: ddot, sum_o, sum_v, norm
!
!     ::::::::::::::::::::::::::::::::::::::::::::::::
!     -::- Construct CNTO transformation matrices -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::::
!
!     Set some local variables
!
      cc2_n_x2am = (wf%n_CC2_o)*(wf%n_CC2_v)*((wf%n_CC2_o)*(wf%n_CC2_v) + 1)/2
      cc2_n_parameters = (wf%n_v)*(wf%n_o) + cc2_n_x2am
      n_cc2_o = wf%n_cc2_o
      n_cc2_v = wf%n_cc2_v
!
!     Open file of CC2 solution vectors
!
      call generate_unit_identifier(unit_solution)
!
!     Determine file name
!
      if (wf%tasks%excited_state) wf%excited_state_specifications%solution_file = 'right_valence'
      if (wf%tasks%core_excited_state) wf%excited_state_specifications%solution_file = 'right_core'
!
      open(unit=unit_solution, file=wf%excited_state_specifications%solution_file,&
            action='read', status='unknown', &
            access='direct', form='unformatted', recl=dp*(cc2_n_parameters), iostat=ioerror) 
!
      if (ioerror .ne. 0) write(unit_output,*) 'Error while opening solution file', ioerror
!
!     -::- Construct M and N -::-
!
      call wf%mem%alloc(M_i_j, n_CC2_o, n_CC2_o)
      call wf%mem%alloc(N_a_b, n_CC2_v, n_CC2_v)
      M_i_j = zero
      N_a_b = zero
!
!     Reading CC2 excitation vectors and summing them
!
      do state = 1, wf%excited_state_specifications%n_singlet_states
!
         call wf%mem%alloc(R, cc2_n_parameters, 1)
!
         read(unit=unit_solution, rec=state) R
!
!        Restict to CC2 active indices (only changes vector if there is a CCS space)
!        This is done so that the CCSD space will always lie inside of the CC2 space
!
         call wf%mem%alloc(R_a_i, n_CC2_v, n_CC2_o)
         call wf%mem%alloc(R_ai_bj, n_CC2_o*n_CC2_v, n_CC2_o*n_CC2_v)
!
         do a = 1, n_CC2_v
            do i = 1, n_CC2_o
!
               ai_full = index_two(a, i, wf%n_v)
!
               R_a_i(a, i) = R(ai_full, 1)
!
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
                     R_ai_bj(ai, bj) = R((wf%n_o)*(wf%n_v) + aibj, 1)
                  enddo
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(R, cc2_n_parameters, 1)
!
!
!        -::- Construct M and N -::-
!
!        Singles contribution
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
                        one,     &
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
                        one,    &
                        N_a_b,   &
                        n_CC2_v)
!
!
         call wf%mem%dealloc(R_a_i, n_cc2_v, n_cc2_o)
!
!        Doubles contribution
!
         call wf%mem%alloc(R_aib_k, (n_CC2_v**2)*(n_CC2_o), n_CC2_o)
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
         call wf%mem%dealloc(R_aib_k, (n_CC2_v**2)*(n_CC2_o), n_CC2_o)
!
         do i = 1, n_CC2_o
            do a = 1, n_CC2_v
               ai = index_two(a, i, n_CC2_v)
               M_i_j(i, i) = M_i_j(i, i) + half*R_ai_bj(ai, ai)
            enddo
         enddo
!
         call wf%mem%alloc(R_a_icj, n_CC2_v, (n_CC2_o**2)*(n_CC2_v))
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
!  
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
         call wf%mem%dealloc(R_a_icj, n_CC2_v, (n_CC2_o**2)*(n_CC2_v))
!
         do a = 1, n_CC2_v
            do i = 1, n_CC2_o
               ai = index_two(a, i, n_CC2_v)
               N_a_b(a, a) = N_a_b(a, a) + half*R_ai_bj(ai, ai)
            enddo
         enddo
!
         call wf%mem%dealloc(R_ai_bj, n_CC2_o*n_CC2_v, n_CC2_o*n_CC2_v)
!
      enddo
!
!     Done with file, delete it
!
      close(unit_solution, status='delete')
!
!     Scale so that M and N are trace 1 matrices. ! OBS OBS OBS denne funker bare med CC2/CCSD tas ut når vi ikke bruker thresholds
!
      call dscal((wf%n_o)*(wf%n_o), one/wf%excited_state_specifications%n_singlet_states, M_i_j, 1)
      call dscal((wf%n_v)*(wf%n_v), one/wf%excited_state_specifications%n_singlet_states, N_a_b, 1)
!
!     -::- Diagonalize M and N matrix -::-
!
      call wf%mem%alloc(eigenvalues_o, n_cc2_o, 1)
      call wf%mem%alloc(work, 4*(n_cc2_o), 1)
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
      write(unit_output,*)'info', info
!
      call wf%mem%dealloc(work, 4*(n_cc2_o), 1)
!
      call wf%mem%alloc(eigenvalues_v, n_cc2_v, 1)
      call wf%mem%alloc(work, 4*(n_cc2_v), 1)
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
      call wf%mem%dealloc(work, 4*(n_cc2_v), 1)
      write(unit_output,*)'M matrix eigenvalues:'
      call vec_print(eigenvalues_o, n_CC2_o, 1)
!
!     -::- Reorder M and N -::-
!
!     dsyev orderes eigenvalues and corresponding eigenvectors in ascending order.
!     We wish to select active space according to highest eigenvalues, thus we must reorder
!
      call wf%mem%alloc(M, n_cc2_o, n_cc2_o)
!
      do i = 1, n_cc2_o
!
         j = i -1
         M(:,i) = M_i_j(:, n_cc2_o - j)
!
      enddo
!
      call wf%mem%dealloc(M_i_j, n_CC2_o, n_CC2_o)
!
      call wf%mem%alloc(N, n_cc2_v, n_cc2_v)
!
      do a = 1, n_cc2_v
!
         b = a -1
         N(:,a) = N_a_b(:,n_cc2_v - b)
!
      enddo
!
      call wf%mem%dealloc(N_a_b, n_CC2_v, n_CC2_v)
!
!     ::::::::::::::::::::::::::::
!     -::- Transform C matrix -::-
!     ::::::::::::::::::::::::::::
!
      call wf%mem%alloc(C_o, wf%n_ao, n_cc2_o)
      call wf%mem%alloc(C_v, wf%n_ao, n_cc2_v)
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
      call wf%mem%alloc(C_o_transformed, wf%n_ao, n_cc2_o)
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
      call wf%mem%dealloc(C_o, wf%n_ao, n_cc2_o)
!
      call wf%mem%alloc(C_v_transformed, wf%n_ao, n_cc2_v)
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
      call wf%mem%dealloc(C_v, wf%n_ao, n_cc2_v)
      call wf%mem%dealloc(N, n_cc2_v, n_cc2_v)
      call wf%mem%dealloc(M, n_cc2_o, n_cc2_o)  
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
      call wf%mem%dealloc(C_o_transformed, wf%n_ao, n_cc2_o)
      call wf%mem%dealloc(C_v_transformed, wf%n_ao, n_cc2_v)
!
!     ::::::::::::::::::::::::::::::::
!     -::- Active space selection -::-
!     ::::::::::::::::::::::::::::::::
!
      sum_o       = 1 - eigenvalues_o(n_CC2_o, 1)
      wf%n_CCSD_o = 1
!
      do while ((sum_o .gt. wf%CCSD_orbitals%delta_o) .and. (wf%n_CCSD_o .lt. n_cc2_o))
!
         sum_o = sum_o - eigenvalues_o(n_CC2_o - (wf%n_CCSD_o), 1)
         wf%n_CCSD_o = wf%n_CCSD_o + 1
!
      enddo
!
      sum_v       = 1 - eigenvalues_v(n_CC2_v, 1)
      wf%n_CCSD_v = 1
!
      do while (sum_v .gt. wf%CCSD_orbitals%delta_v .and. (wf%n_CCSD_v .lt. n_cc2_v))
!
         sum_v = sum_v - eigenvalues_v(n_CC2_v - (wf%n_CCSD_v), 1)
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
      call wf%mem%dealloc(eigenvalues_o, n_cc2_o, 1)
      call wf%mem%dealloc(eigenvalues_v, n_cc2_v, 1)
!
!     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     -::- Finding orbital energies and new block diagonal C matrix -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     Initialize amplitudes and associated attributes
!
      call wf%initialize_single_amplitudes
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf%read_transform_cholesky
!
!     Initialize fock matrix
!
      call wf%initialize_fock_matrix
!
!     -::- Occupied orbitals -::-
!
!     Diagonalize active-active block
!
      if (wf%n_CCSD_o .gt. 0) then
!
         call wf%mem%alloc(work, 4*(wf%n_CCSD_o), 1)
         call wf%mem%alloc(orbital_energies, (wf%n_CCSD_o), 1)
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
         call wf%mem%dealloc(work, 4*(wf%n_CCSD_o), 1)
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
         call wf%mem%dealloc(orbital_energies, (wf%n_CCSD_o), 1)
      endif
!
!     Diagonalize inactive-inactive block 
!
      if (wf%n_CC2_o .gt. 0) then
         call wf%mem%alloc(work, 4*wf%n_CC2_o, 1)
         call wf%mem%alloc(orbital_energies, wf%n_CC2_o, 1)
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
         call wf%mem%dealloc(work, 4*wf%n_CC2_o, 1)
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
         call wf%mem%dealloc(orbital_energies, wf%n_CC2_o, 1)
      endif
!
!     Transform C-matrix to block diagonal (active occupied block)
!
      call wf%mem%alloc(C_o, wf%n_ao, wf%n_o)
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
         call wf%mem%alloc(C_o_transformed, wf%n_ao, wf%n_CCSD_o)
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
         call wf%mem%dealloc(C_o_transformed, wf%n_ao, wf%n_CCSD_o)
      endif
!
!     Transform C-matrix to block diagonal (inactive occupied block)
!
      if (wf%n_CC2_o .gt. 0) then
         call wf%mem%alloc(C_o_transformed, wf%n_ao, wf%n_CC2_o)
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
         call wf%mem%dealloc(C_o, wf%n_ao, wf%n_o)
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
         call wf%mem%dealloc(C_o_transformed, wf%n_ao, wf%n_CC2_o)
      endif
!
!     -::- Vacant orbitals -::-
!
!     Diagonalize active-active block
!
      if (wf%n_CCSD_v .gt. 0) then
         call wf%mem%alloc(work, 4*(wf%n_CCSD_v), 1)
         call wf%mem%alloc(orbital_energies, (wf%n_CCSD_v), 1)
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
         call wf%mem%dealloc(work, 4*(wf%n_CCSD_v), 1)
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
         call wf%mem%dealloc(orbital_energies, (wf%n_CCSD_v), 1)
!
      endif
!
!     Diagonalize inactive-inactive block 
!
      if (wf%n_CC2_v .gt. 0) then
         call wf%mem%alloc(work, 4*wf%n_CC2_v, 1)
         call wf%mem%alloc(orbital_energies, wf%n_CC2_v, 1)
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
         call wf%mem%dealloc(work, 4*wf%n_CC2_v, 1)
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

         call wf%mem%dealloc(orbital_energies, wf%n_CC2_v, 1)
!
      endif
!
!     Transform C-matrix to block diagonal (active virtual block)
!
      call wf%mem%alloc(C_v, wf%n_ao, wf%n_v)
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
         call wf%mem%alloc(C_v_transformed, wf%n_ao, wf%n_CCSD_v)
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
!
               ij = index_two(i, j + wf%n_o, wf%n_ao)
               wf%mo_coef(ij, 1) = C_v_transformed(i, j) 
!
            enddo
         enddo
!
         call wf%mem%dealloc(C_v_transformed, wf%n_ao, wf%n_CCSD_v)
!
      endif
!
!     Transform C-matrix block diagonal (inactive virtual block)
!
      if (wf%n_CC2_v .gt. 0) then
!
         call wf%mem%alloc(C_v_transformed, wf%n_ao, wf%n_CC2_v)
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
         call wf%mem%dealloc(C_v_transformed, wf%n_ao, wf%n_CC2_v)
!
      endif
!
      call wf%mem%dealloc(C_v, wf%n_ao, wf%n_v)     
!
      call wf%destruct_single_amplitudes
!
   end subroutine cnto_orbitals_mlccsd
!
!
   module subroutine cnto_init_cc2_mlccsd(wf, cc2_wf)
!!
!!    Initialize lower level method,
!!    Written by Sarai D. Folkestad, Dec 2017  
!!
!!    Initializes 
!!
      implicit none
!
      class(mlccsd) :: wf
!
      class(mlcc2)    :: cc2_wf
!
      integer(i15) :: lower_level_n_singlet_states
      integer(i15) :: i = 0, j = 0
!
      integer(i15), dimension(:,:), allocatable :: index_list
!
      logical :: start_vec_exists
!
!     Set calculation tasks
!
      cc2_wf%tasks = wf%tasks
!
!     Set calculation settings
!
      cc2_wf%settings = wf%settings
!
      call cc2_wf%mem%init(int(wf%mem%available/1.0D9))
!
      cc2_wf%core_excited_state_specifications  = wf%core_excited_state_specifications
      cc2_wf%excited_state_specifications       = wf%excited_state_specifications
!
!     Set convergence threshold for lower lying method
!
      cc2_wf%ground_state_specifications%energy_threshold = 1.0D-08 
      cc2_wf%ground_state_specifications%residual_threshold = 1.0D-08 
!
      cc2_wf%excited_state_specifications%energy_threshold = 1.0D-04 
      cc2_wf%excited_state_specifications%residual_threshold = 1.0D-04  
!
!     Set mlcc settings
!
      cc2_wf%mlcc_settings = wf%mlcc_settings
      cc2_wf%mlcc_settings%CCSD = .false.
!
      cc2_wf%CC2_orbitals = wf%CC2_orbitals
!
!     Test for user specified start vector
!
      if (wf%excited_state_specifications%user_specified_start_vector) then
!
!        Since orbitals will swap order, start vector in higher level method must be removed
!
         wf%excited_state_specifications%user_specified_start_vector = .false.
         call wf%mem%dealloc_int(wf%excited_state_specifications%start_vectors, wf%excited_state_specifications%n_singlet_states, 1)
!        
      endif
!
!     :: Initialize lower level method ::
!  
      cc2_wf%name = 'MLCC2'
!
!     Set implemented generic methods
!
      cc2_wf%implemented%ground_state         = .true.
      cc2_wf%implemented%excited_state        = .true.
      cc2_wf%implemented%core_excited_state   = .true.
!
!     Read Hartree-Fock info from SIRIUS
!
      call cc2_wf%read_hf_info
!
      if (cc2_wf%mlcc_settings%CCS) then
!
         call cc2_wf%orbital_partitioning
!
      else
!
!        Do full space CC2 calculation
!
         cc2_wf%n_CC2_o = wf%n_o
         cc2_wf%n_CC2_v = wf%n_v
!
         cc2_wf%first_CC2_o = 1
         cc2_wf%first_CC2_v = 1         
!
      endif
!
!     Initialize amplitudes and associated attributes
!
      call cc2_wf%initialize_single_amplitudes
!
!     Set the number of parameters in the wavefunction
!     (that are solved for in the ground and excited state solvers) 
!
      cc2_wf%n_t1am = (cc2_wf%n_o)*(cc2_wf%n_v)
      cc2_wf%n_parameters = cc2_wf%n_t1am
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call cc2_wf%read_transform_cholesky
!
!     Initialize fock matrix
!
      call cc2_wf%initialize_fock_matrix
      call cc2_wf%destruct_single_amplitudes
!
   end subroutine cnto_init_cc2_mlccsd
!
!
   module subroutine print_orbital_info_mlccsd(wf)
!!
!!    Print CNTO info, 
!!    Written by Sarai D. Folkestad, Aug. 2017
!!
!!    Prints information on CNTO partitioning
!!
      implicit none 
!
      class(mlccsd) :: wf
!
      write(unit_output, '(t6,a11)')'CCSD space:'
      write(unit_output,'(t6,a34, i3)') 'Number of occupied orbitals:      ', wf%n_CCSD_o
      write(unit_output,'(t6,a34, i3/)')'Number of virtual orbitals:       ', wf%n_CCSD_v
!
      write(unit_output, '(t6,a10)')'CC2 space:'
      write(unit_output,'(t6,a34, i3)') 'Number of occupied orbitals:      ', wf%n_CC2_o
      write(unit_output,'(t6,a34, i3/)')'Number of virtual orbitals:       ', wf%n_CC2_v
!  
      write(unit_output, '(t6,a10)')'CCS space:'
      write(unit_output,'(t6,a34, i3)') 'Number of occupied orbitals:      ', wf%n_CCS_o
      write(unit_output,'(t6,a34, i3/)')'Number of virtual orbitals:       ', wf%n_CCS_v
      flush(unit_output)
!
   end subroutine print_orbital_info_mlccsd
!
end submodule orbital_partitioning
