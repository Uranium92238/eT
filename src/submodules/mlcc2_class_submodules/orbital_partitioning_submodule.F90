submodule (mlcc2_class) orbital_partitioning
!
!!
!!    Orbital partitioning submodule (MLCC2) 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2017
!!
!!    Contains the following family of procedures of the MLCC2 class:
!!
!!    orbital_partitioning:         Directs the orbital partitioning
!!    cholesky_localization_drv:    Directs orbital localization by cholesky decomposition  
!!    cholesky_orbital_constructor: Directs construction of new orbitals                     
!!    cholesky_decomposition:       Cholesky decomposes the density
!!    cholesky_orbitals:            Constructs new orbitals (C matrix) from cholesky vectors 
!! 
!!    Contains the following module subroutines and functions: These should eventually be moved to some 
!!                                                       utils !
!!
!!    get_number_of_active_spaces:     Returns number of active spaces from cholesky.inp
!!    get_number_of_active_atoms:      Returns number of active atoms from cholesky.inp
!!    get_active_atoms:                Returns the indices of the active atoms from cholesky.inp
!!    construct_active_ao_index_list:  Returns active_ao_index_list to be used as pivoting elements in cholesky decomposition
!!    read_atom_info:                  Reads info from dalton, n_nuclei and n_ao
!!    read_center_info:                Reads info from dalton, which aos belong to which nuclei
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
   module subroutine orbital_partitioning_mlcc2(wf)
!!
!!    Orbital partitioning,
!!    Written by Sarai D. Folkestad, June 2017
!!
!!    Directs the partitioning for mlcc calculations.
!!
!!    So far only Cholesky decomposition is available. 
!!
      implicit none
!
      class(mlcc2) :: wf
!
      if (wf%mlcc_settings%cholesky) then
!
!        If cholesky - do Cholesky decomposition
!
         call wf%cholesky_localization_drv
!
      elseif (wf%mlcc_settings%cnto) then
!
         call wf%cnto_orbital_drv
         call wf%print_cnto_info
!
      endif
!
   end subroutine orbital_partitioning_mlcc2
!
!
   module subroutine cholesky_localization_drv_mlcc2(wf)
!!
!!    Cholesky orbital localization. driver,
!!    Written by Sarai D. Folkestad, June 2017
!!
!!    Driver for Cholesky density decomposition.  
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
      class(mlcc2) :: wf
!
      real(dp), dimension(:,:), allocatable :: orbitals
      real(dp), dimension(:,:), allocatable :: orbital_energies
!  
      integer(i15) :: n_nuclei
      integer(i15) :: n_CC2_atoms = 0
      integer(i15) :: offset_o = 0, offset_v = 0
      integer(i15) :: n_active_aos = 0
!
      integer(i15), dimension(:,:), allocatable :: ao_center_info, n_ao_on_center
      integer(i15), dimension(:,:), allocatable :: active_atoms
      integer(i15), dimension(:,:), allocatable :: active_ao_index_list
      integer(i15) :: n_vectors_o = 0, n_vectors_v = 0
!
!     Timing variables
!
      real(dp) :: start_chol_deco = 0, end_chol_deco = 0
!
!     IO-variables
!
      logical      :: file_exists     
      integer(i15) :: unit_cholesky_decomp = 0, unit_overlap = 0, ioerror = 0
!
!     Indices
!
      integer(i15) :: i = 0, j = 0, ij = 0
!
!     Allocatables
!
      real(dp), dimension(:,:), allocatable :: density_o, density_v
      real(dp), dimension(:,:), allocatable :: ao_fock
      real(dp), dimension(:,:), allocatable :: C
      real(dp), dimension(:,:), allocatable :: S_packed
      real(dp), dimension(:,:), allocatable :: S
      real(dp), dimension(:,:), allocatable :: Y_1, Y_2
!
!     Start timings
!
      call cpu_time(start_chol_deco)
!
!     Prepare for decomposition
!
!     Prints
!
      write(unit_output,'(/t3,a/)')    ':: Cholesky decomposition '
!
!     :: Get center info ::
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
!     :::::::::::::::::::::::
!     -::- Occupied part -::-
!     :::::::::::::::::::::::
!
!     :: Construct AO fock matrix ::
!
      call allocator(ao_fock, wf%n_ao, wf%n_ao)
      call wf%construct_ao_fock_new(ao_fock)
!
      call allocator(density_o, wf%n_ao, wf%n_ao)
      density_o = zero     
!
      call wf%construct_density_matrix(density_o, wf%mo_coef, wf%n_o, wf%n_v)
!
!     Prepare for loop over active spaces
!
      offset_o = 1
!
!     Variables for storing information on spaces
!     Allocations for Cholesky localized orbitals
!
      call allocator (orbitals, wf%n_ao, wf%n_mo)
      call allocator (orbital_energies, wf%n_mo, 1)
!
!     Get CC2-active atoms
!
      n_CC2_atoms =  get_number_of_active_atoms(unit_cholesky_decomp, 'CC2  ')
!
!     Allocate active atoms list
!
      call allocator_int(active_atoms, n_CC2_atoms, 1)
!
!     Get list of active atoms
!
      call get_active_atoms(unit_cholesky_decomp, active_atoms,  n_CC2_atoms, 'CC2  ')
!
!     Sanity check on active atoms
!
      if (n_CC2_atoms .gt. n_nuclei) then
         write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
         stop
      endif
!
      do i = 1, n_CC2_atoms
         if (active_atoms(i,1) .gt. n_nuclei) then
            write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
            stop
         endif
      enddo
!
!     :: Constructing active (CC2) localized Cholesky orbitals ::
!
      n_active_aos = 0
      do i = 1, n_CC2_atoms
         n_active_aos = n_active_aos + n_ao_on_center(active_atoms(i,1),1)
      enddo
!
!     Construct active_ao_index_list
!
      call allocator_int(active_ao_index_list, n_active_aos, 1)
      call construct_active_ao_index_list(active_ao_index_list, n_active_aos, active_atoms, &
                                               n_CC2_atoms, ao_center_info, wf%n_ao)
!
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                           density_o, n_vectors_o,&
                                          .true., n_active_aos, active_ao_index_list)
!
!     Save active space information
!
      wf%n_CC2_o = n_vectors_o
!
      wf%first_CC2_o = offset_o
!
!     Calculate new offset         
!
      offset_o = offset_o + n_vectors_o
!
      call deallocator_int(active_atoms, n_CC2_atoms, 1)
      call deallocator_int(active_ao_index_list, n_active_aos, 1)
!
!     :: CCS  localized Cholesky orbitals  ::
!
      n_vectors_o = 0
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_o, ao_fock, &
                                           density_o, n_vectors_o,&
                                           .false., n_active_aos)
      wf%n_CCS_o     = n_vectors_o
!
      wf%first_CCS_o = offset_o
!
      call deallocator(density_o, wf%n_ao, wf%n_ao)
!
!     ::::::::::::::::::::::
!     -::- Virtual part -::-
!     ::::::::::::::::::::::

!
!     Allocations for Cholesky localized orbitals
!   
      call allocator(density_v, wf%n_ao, wf%n_ao)
      density_v = zero
!
      call wf%construct_density_matrix_v(density_v, wf%mo_coef, wf%n_o, wf%n_v)
!
      offset_v = 1 + wf%n_o
!
!     Get CC2-active atoms
!
      n_CC2_atoms =  get_number_of_active_atoms(unit_cholesky_decomp, 'CC2  ')
!
!     Allocate active atoms list
!
      call allocator_int(active_atoms, n_CC2_atoms, 1)
!
!     Get list of active atoms
!
      call get_active_atoms(unit_cholesky_decomp, active_atoms,  n_CC2_atoms, 'CC2  ')
!
!     Sanity check on active atoms
!
      if (n_CC2_atoms .gt. n_nuclei) then
         write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
         stop
      endif
!
      do i = 1, n_CC2_atoms
         if (active_atoms(i,1) .gt. n_nuclei) then
            write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
            stop
         endif
      enddo
!
!     :: Constructing active (CC2) localized Cholesky orbitals ::
!
      n_active_aos = 0
      do i = 1, n_CC2_atoms
         n_active_aos = n_active_aos + n_ao_on_center(active_atoms(i,1),1)
      enddo
!
!     Construct active_ao_index_list
!
      call allocator_int(active_ao_index_list, n_active_aos, 1)
      call construct_active_ao_index_list(active_ao_index_list, n_active_aos, active_atoms, &
                                                  n_CC2_atoms, ao_center_info, wf%n_ao)
!
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                              density_v, n_vectors_v, &
                                              .true., n_active_aos, active_ao_index_list)
!
!     Save active space information
!
      wf%n_CC2_v = n_vectors_v
!
      wf%first_CC2_v = offset_v - wf%n_o
!
!     Calculate new offset         
!
      offset_v = offset_v + n_vectors_v
!
      call deallocator_int(active_atoms, n_CC2_atoms, 1)
      call deallocator_int(active_ao_index_list, n_active_aos, 1)

!
!
!     :: CCS ::
!
      n_vectors_v = 0
      call wf%cholesky_orbital_constructor(orbitals, orbital_energies, offset_v, ao_fock, &
                                           density_v, n_vectors_v,&
                                           .false., n_active_aos)
!
      call deallocator(density_v, wf%n_ao, wf%n_ao)     
      call deallocator(ao_fock, wf%n_ao, wf%n_ao)
!
!     Save inactive space information
!
      wf%n_CCS_v = n_vectors_v
!
      wf%first_CCS_v = offset_v - wf%n_o
!
      wf%mo_coef        = orbitals
      wf%fock_diagonal  = orbital_energies 
!
      call deallocator (orbitals, wf%n_ao, wf%n_mo)
      call deallocator (orbital_energies, wf%n_mo, 1)  
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
!     :: Save new orbitals and orbital energies ::
!
!
      call deallocator_int(n_ao_on_center, n_nuclei, 2)
      call deallocator_int(ao_center_info, wf%n_ao, 2)
!     
!     Close cholesky.inp
!
      close(unit_cholesky_decomp)
!
!     Print decomposition info
!
      write(unit_output, '(/t3,a15, i3/)')'Active space:  '
      write(unit_output,'(t3,a40, i3)') 'Number of active occupied orbitals:  ', wf%n_CC2_o
      write(unit_output,'(t3,a40, i3/)')'Number of active virtual orbitals:   ', wf%n_CC2_v
!
      write(unit_output, '(/t3,a15/)')'Inactive space:  '
      write(unit_output,'(t3,a40, i3)') 'Number of inactive occupied orbitals:  ', wf%n_CCS_o
      write(unit_output,'(t3,a40, i3/)')'Number of inactive virtual orbitals:   ', wf%n_CCS_v
!
!     Print timings
!
      call cpu_time(end_chol_deco)
      if (timings) write(unit_output,'(/t3,a50,f14.8/)') 'Total time in Cholesky orbital construction (seconds):',&
                            end_chol_deco - start_chol_deco
      flush(unit_output)
!
!
   end subroutine cholesky_localization_drv_mlcc2
!
!
   module subroutine cholesky_orbital_constructor_mlcc2(wf, orbitals, orbital_energies, offset, ao_fock, density, n_vectors,&
                              selection, n_active_aos, active_ao_index_list)
!!
!!    Cholesky orbital constructor,
!!    Written by Sarai Dery Folkestad, June 2017
!!
!!    Constructs new localized orbitals (occupied/virtual) by 
!!    - Decomposing the density (occupied/virual)
!!    - Transforming Fock matrix with Cholesky vectors, and diagonalizing it.
!!      New orbitals are eigenvectors, orbital energies are eigenvectors.
!!
      implicit none
!
      class(mlcc2)                                       :: wf
      real(dp), dimension(wf%n_ao, wf%n_mo)              :: orbitals
      real(dp), dimension(wf%n_mo, 1)                    :: orbital_energies
      real(dp), dimension(wf%n_ao, wf%n_ao)              :: ao_fock
      integer(i15)                                       :: n_active_aos, offset
      integer(i15)                                       :: n_vectors
      real(dp), dimension(wf%n_ao,wf%n_ao)               :: density
      logical                                            :: selection
      integer(i15), dimension( n_active_aos,1), optional :: active_ao_index_list
!
      real(dp), dimension(:,:), allocatable              :: cholesky
      real(dp), dimension(:,:), allocatable              :: orbitals_temp
      real(dp), dimension(:,:), allocatable              :: orbital_energies_temp
!
      integer(i15) :: i = 0, j = 0
!
      if (selection .and. (.not. present(active_ao_index_list)) ) then
         write(unit_output,*) 'WARNING: Illegal argument list for cholesky decomposition'
            stop
      endif
!
      call allocator(cholesky, wf%n_ao, wf%n_ao)     
        
      cholesky = zero   
!
!     :: Construct cholesky vectors by decomposing density ::
!
      if (selection) then
         call wf%cholesky_decomposition(density, cholesky, &
                                               n_vectors, .true., n_active_aos, active_ao_index_list)
      else 
         call wf%cholesky_decomposition(density, cholesky, &
                                               n_vectors, .false., n_active_aos)
      endif   
!
!     :: Construct cholesky orbitals ::
!
      call allocator(orbitals_temp, wf%n_ao, n_vectors)
      call allocator(orbital_energies_temp, n_vectors,1)
      orbitals_temp = zero     
      orbital_energies_temp = zero
!
      call wf%cholesky_orbitals(cholesky, n_vectors, orbitals_temp, orbital_energies_temp, ao_fock)       
!
!     Place new orbitals and orbital energies in array
!
      do i = 1, wf%n_ao
         do j = 1, n_vectors
            orbitals(i,j + offset - 1)         = orbitals_temp(i, j)
            orbital_energies(j + offset - 1,1) = orbital_energies_temp(j,1)
         enddo
      enddo
!
!     Deallocations
!
      call deallocator(cholesky, wf%n_ao, wf%n_ao)
      call deallocator(orbitals_temp, wf%n_ao, n_vectors)
      call deallocator(orbital_energies_temp, n_vectors,1)
!
   end subroutine cholesky_orbital_constructor_mlcc2
!
!
   module subroutine cholesky_decomposition_mlcc2(wf, density, cholesky_vectors,&
                                                     n_vectors, selection, n_active_aos, active_ao_index_list)
!!
!!    Cholesky decomposition, 
!!    Written by Sarai dery Folkestad, June 2017.
!!
!!    Cholesky decomposes the density (occupied/virtual).
!!    Pivoting elements are chosen according to the active ao-list if it is pressent.
!!    If not, maximum diagonal elements are chosen as pivoting elements.
!!
!!    The Cholesky vectors are subtracted from the incoming density matrix, and the returned density can be furteh decomposed
!!    for inactive region.
!! 
      implicit none
!
      class(mlcc2)                                       :: wf
      integer(i15)                                       :: n_active_aos
      integer(i15)                                       :: n_vectors
      real(dp), dimension(wf%n_ao,wf%n_ao)               :: density
      real(dp), dimension(wf%n_ao, wf%n_ao)              :: cholesky_vectors
      logical                                            :: selection
      integer(i15), dimension( n_active_aos,1), optional :: active_ao_index_list
!
      integer(i15) :: index_max = 0
      real(dp)     :: max_diagonal
!
      real(dp), parameter  :: threshold = 1.0d-02
      real(dp), parameter  :: tolerance = 1.0d-09
!
      integer(i15) :: i = 0, j = 0, k = 0
!
      if (selection .and. (.not. present(active_ao_index_list)) ) then
         write(unit_output,*) 'WARNING: Illegal argument list for cholesky decomposition'
            stop
      endif
!
!     Looping over the number of cholesky vectors   
!
      do i = 1, wf%n_ao
         n_vectors = i
!
!        Find the maximum diagonal
!
         index_max = 0
         max_diagonal = -1.0d10 ! Make this a really small number
!
         if (selection) then
            do j = 1, n_active_aos
!
               if (density(active_ao_index_list(j,1),active_ao_index_list(j,1)) .gt. max_diagonal) then
!
                  max_diagonal = density(active_ao_index_list(j,1),active_ao_index_list(j,1))
                  index_max    = active_ao_index_list(j,1)
!
               endif
!
            enddo
         else
             do j = 1, wf%n_ao
!
               if (density(j,j) .gt. max_diagonal) then
!
                  max_diagonal = density(j,j)
                  index_max    = j
!
               endif
!
            enddo
         endif
!
!        Check against threshold and whether diagonal is negative
!
         if (max_diagonal .lt. 0.0d0) then
            if (abs(max_diagonal) .gt. tolerance) then
!
               write(unit_output,*)'WARNING: Found negative diagonal in cholesky decomposition.'
               stop
!
            endif
         endif
!
         if (abs(max_diagonal) .lt. threshold) then
            n_vectors = n_vectors - 1                
            return
         endif
!
!        Cholesky vectors
!
         do j = 1, wf%n_ao
!
            cholesky_vectors(j,i) = density(j, index_max)/sqrt(max_diagonal)
!
         enddo
!
!        Subtract from density
!
         do j = 1, wf%n_ao
            do k = 1, wf%n_ao
!
               density(k,j) = density(k,j) - cholesky_vectors(k,i)*cholesky_vectors(j,i)
!
            enddo
         enddo
!
         do j = 1, wf%n_ao
            density(j,index_max) = 0.0D0
            density(index_max,j) = 0.0D0
         enddo
!
      enddo

   end subroutine cholesky_decomposition_mlcc2
!
!
   module subroutine cholesky_orbitals_mlcc2(wf, cholesky_vectors, n_vectors, orbitals, orbital_energies, ao_fock)
!!
!!    Cholesky orbitals,
!!    Written by Sarai Dery Folkestad, June 2017
!!
!!    Makes the new MO's fromthe Cholesky vectors   
!!    - Transforms the AO fock matrix by the Cholesky vectors
!!    - Diagonalize the MO fock to get orbital energies and new orbitals.
!!
      implicit none
!
      class(mlcc2)                              :: wf
      real(dp), dimension(wf%n_ao, wf%n_ao)     :: cholesky_vectors
      real(dp), dimension(wf%n_ao, n_vectors)   :: orbitals
      real(dp), dimension(n_vectors, 1)         :: orbital_energies
      integer(i15)                              :: n_vectors
      real(dp), dimension(wf%n_ao,wf%n_ao)      :: ao_fock
!
      real(dp), dimension(:,:), allocatable :: mo_fock
      real(dp), dimension(:,:), allocatable :: X
      real(dp), dimension(:,:), allocatable :: work
      integer(i15) :: lwork
      integer(i15) :: info = 0
!
      integer(i15) :: a = 0, b = 0, i = 0, j = 0
!
     if (n_vectors == 0) return
!
      call allocator(X, n_vectors, wf%n_ao)
!
!     :: Transform AO Fock :: 
!
      call dgemm('T', 'N',          &
                  n_vectors,        &
                  wf%n_ao,          &
                  wf%n_ao,          &
                  one,              &
                  cholesky_vectors, &
                  wf%n_ao,          &
                  ao_fock,          &
                  wf%n_ao,          &
                  zero,             &
                  X,                &
                  n_vectors)
!
      call allocator(mo_fock, n_vectors, n_vectors)
!
      call dgemm('N', 'N',          &
                  n_vectors,        &
                  n_vectors,        &
                  wf%n_ao,          &
                  one,              &
                  X,                &
                  n_vectors,        &
                  cholesky_vectors, &
                  wf%n_ao,          &
                  zero,             &
                  mo_fock,          &
                  n_vectors)
!
      call deallocator(X, n_vectors, wf%n_ao)
!
!     :: Diagonalize MO fock ::
!
      call allocator(work, 4*n_vectors, 1)
      work = zero
!
      call dsyev('V','U', &
                  n_vectors, &
                  mo_fock, &
                  n_vectors, &
                  orbital_energies, &
                  work, & 
                  4*n_vectors, &
                  info)
!
      call deallocator(work, 4*n_vectors, 1)
!
      if (info .ne. 0) then
         write(unit_output,*)'WARNING: Diagonalization not successful.'
               stop
      endif
!
!     Construct C-matrix (orbital coefficients)
!
      call dgemm('N', 'N',          &
                  wf%n_ao,          &  
                  n_vectors,        &
                  n_vectors,        &
                  one,              &
                  cholesky_vectors, &
                  wf%n_ao,          &
                  mo_fock,          &
                  n_vectors,        &
                  zero,             &
                  orbitals,         &
                  wf%n_ao)
!
!
      call deallocator(mo_fock, n_vectors, n_vectors)
!
   end subroutine cholesky_orbitals_mlcc2
!
!
   function get_number_of_active_atoms(unit_cholesky_decomp, ml_level)
!!
!!    Get number of active atoms
!!    Written by Sarai D. Folkestad June 2017
!!
!!    Reads cholesky.inp, and returns the number of atoms treated at ml_level of the CC hierarchy
!!    in active space of question.
!!
      implicit none
!
      integer(i15)      :: get_number_of_active_atoms
      integer(i15)      :: unit_cholesky_decomp 
      character(len=5)  :: ml_level             ! CC2/CCSD/CC3
!
      character(len=40) :: line
      integer(i15)      :: ioerror = 0
      integer(i15)      :: previous_spaces, previous_atoms, previous_n_atoms
!
      rewind(unit_cholesky_decomp)
!
      do
         read(unit_cholesky_decomp,'(a40)')line
!
         do while (line(1:1) == '!' .or. trim(line) == '') ! Comment or blank line: read the next line
!
            read(unit_cholesky_decomp,'(a40)') line
         enddo
!  
         if (trim(line) == ml_level) then
!
            read(unit_cholesky_decomp,'(i5)', iostat=ioerror) get_number_of_active_atoms
!
            if (ioerror .ne. 0) then
               write(unit_output,*)'WARNING: Error while reading number of active spaces from cholesky.inp'
               stop
            endif
            exit
!
         elseif(trim(line) == '#end of Cholesky input') then
!
            backspace(unit_cholesky_decomp)
            get_number_of_active_atoms = 0
            return
!
         endif
      enddo
!
   end function get_number_of_active_atoms
!
!
   module subroutine get_active_atoms(unit_cholesky_decomp, active_atoms, n_active_atoms ,ml_level)
!!
!!    Get active atoms
!!    Written by Sarai D. Folkestad June 2017
!!
!!    Reads cholesky.inp, and returns the indices of the active atoms treated at ml_level of the CC hierarchy
!!    in active space of question.
!!
      implicit none
!
      integer(i15)      :: unit_cholesky_decomp
      integer(i15)      :: n_active_atoms
      character(len=5)  :: ml_level
!
      integer(i15), dimension(n_active_atoms,1) :: active_atoms
!
      character(len=40) :: line
      integer(i15)      :: ioerror = 0
      integer(i15)      :: n_atoms
!
      rewind(unit_cholesky_decomp)
      do
         read(unit_cholesky_decomp,'(a40)')line
!
         do while (line(1:1) == '!' .or. trim(line) == '') 
!
            read(unit_cholesky_decomp,'(a40)') line
         enddo
!  
         if (trim(line) == ml_level) then
!
            read(unit_cholesky_decomp,*, iostat=ioerror) n_atoms
!
            read(unit_cholesky_decomp,*, iostat=ioerror) active_atoms
!
            if (ioerror .ne. 0) then
               write(unit_output,*)ioerror
               write(unit_output,*)'WARNING: Error while reading number of active atoms from cholesky.inp'
               stop
            endif
            exit
!
         elseif(trim(line) == '#end of Cholesky input') then
!
            backspace(unit_cholesky_decomp)
            return
!
         endif
      enddo
   end subroutine get_active_atoms
!
!
      module subroutine construct_active_ao_index_list(active_ao_index_list, n_active_aos, active_atoms, &
                                                   n_active_atoms, ao_center_info, n_ao)
!!
!!    Construct active ao index list,
!!    Written by Sarai Dery Folkestad, June 2017.
!!
!!    Constructs list of active ao's for cholesky decomposition.
!!
      implicit none
!
      integer(i15)                               :: n_active_aos
      integer(i15)                               :: n_active_atoms
      integer(i15)                               :: n_ao
      integer(i15), dimension(n_active_aos, 1)   :: active_ao_index_list
      integer(i15), dimension(n_active_atoms, 1) :: active_atoms
      integer(i15), dimension(n_ao, 2)           :: ao_center_info
!
      integer(i15) :: i = 0, j = 0, counter = 0
!
      counter = 0    
      do i = 1, n_active_atoms 
         do j = 1, n_ao
            if (ao_center_info(j,1) == active_atoms(i,1)) then
               counter = counter + 1
               active_ao_index_list(counter,1) = ao_center_info(j,2)
            endif
         enddo
      enddo
!
    end subroutine construct_active_ao_index_list
!
!
   module subroutine read_atom_info(n_nuclei, n_ao)
!!
!!    Read atom info,
!!    Written by Sarai Dery Folkestad, June 2017.
!!
!!    Reads atom info from DALTON generated file:
!!    Reads:
!!       - Number of nuclei
!!       - Number of AO's
!!
      implicit none
!
      integer(i15) :: n_nuclei,n_ao 
!
      integer(i15) :: unit_center = 0
      integer(i15) :: ioerror     = 0
!
      call generate_unit_identifier(unit_center)
      open(unit=unit_center, file='center_info', status='unknown', form='unformatted', iostat=ioerror)
      if (ioerror .ne. 0) write(unit_output,*)'WARNING: Error while opening center_info'
      rewind(unit_center)
!
!     Read number of nuclei and aos
!
      read(unit_center) n_nuclei, n_ao
!
      close(unit_center)
!
   end subroutine read_atom_info
!
   module subroutine read_center_info(n_nuclei, n_ao, n_ao_on_center, ao_center_info)
!!
!!    Read center info,
!!    Written by Sarai Dery Folkestad, June 2017.
!!
!!    Reads atom info from DALTON generated file:
!!    Reads:
!!       - Information of which ao's belong to which nuclei
!!
      implicit none
!
      integer(i15) :: n_nuclei
      integer(i15) :: n_ao
      integer, dimension(n_nuclei, 1)  :: n_ao_on_center
      integer, dimension(n_ao, 2)      :: ao_center_info
!
      integer(i15) :: offset = 0, nucleus = 0
      integer(i15) :: i = 0
      integer(i15) :: unit_center = 0, ioerror = 0
!
!     Read number of aos on each center
!
      call generate_unit_identifier(unit_center)
      open(unit=unit_center, file='center_info', status='unknown', form='unformatted', iostat=ioerror)
      if (ioerror .ne. 0) write(unit_output,*)'WARNING: Error while opening center_info'
      rewind(unit_center)
!
!     Empty read 
!
      read(unit_center)
!
      read(unit_center, iostat=ioerror) n_ao_on_center
      if (ioerror .ne. 0) write(unit_output,*)'WARNING: Error while reading center_info'
!
      offset = 1
!
      do nucleus = 1, n_nuclei
!
         do i = 1, n_ao_on_center(nucleus, 1)
            ao_center_info(offset + i - 1,1)=nucleus
         enddo
!
!        Read which aos are on which centers
!
         read(unit_center, iostat=ioerror) (ao_center_info(offset + i - 1, 2), i = 1, n_ao_on_center(nucleus, 1))
         if (ioerror .ne. 0) write(unit_output,*)'WARNING: Error while reading center_info'
!
         offset = offset + n_ao_on_center(nucleus,1)
!
      enddo
!
      close(unit_center)
!
   end subroutine read_center_info
!
!     :::::::::::::::::::::::
!     -::- CNTO Routines -::-
!     :::::::::::::::::::::::
!
   module subroutine cnto_orbital_drv_mlcc2(wf)
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
      class(mlcc2) :: wf
!
      real(dp) :: start_cnto = 0, end_cnto = 0
!     
!     Prints
!     
      write(unit_output,'(/t3,a/)') ':: CNTO orbital partitioning for MLCC2 calculation '
!
!     Timings 
!
      call cpu_time(start_cnto)
!
!     CNTO orbital selection for CC2
!
      call wf%cc2_cnto
!
!     Timings
!
      call cpu_time(end_cnto)
!
      if (timings) write(unit_output,'(/t3,a50,f14.8/)') 'Total time in cnto orbital construction (seconds):', end_cnto - start_cnto
      flush(unit_output)
!
   end subroutine cnto_orbital_drv_mlcc2
!
!
   module subroutine cc2_cnto_mlcc2(wf)
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
      class(mlcc2) :: wf
!
      type(ccs), allocatable :: ccs_wf
!
      integer(i15) :: unit_solution = -1
      integer(i15) :: ioerror = 0
      integer(i15) :: lower_level_n_singlet_states
      integer(i15) :: state
!
      real(dp), dimension(:,:), allocatable :: solution
      real(dp), dimension(:,:), allocatable :: R_a_i
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
      integer(i15) :: i = 0, j = 0, a = 0, b = 0, ai = 0, ij = 0
!
      real(dp) :: trace, ddot, sum_o, sum_v
!
!     ::::::::::::::::::::::::::::::::::::::::::::::::
!     -::- Running lower level method calculation -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::::
!
!     Allocate lower level method
!
      allocate(ccs_wf)
!
!     Set calculation tasks
!
      ccs_wf%tasks%ground_state = .true.
      ccs_wf%tasks%excited_state = .true.
!
!     Set calculation settings
!
      ccs_wf%settings = wf%settings
! 
!     Set number of excitations to use for cnto generation
!        - Should be some function of wf%tasks%n_singlet_states@
!
      lower_level_n_singlet_states = 4*wf%tasks%n_singlet_states
! 
      ccs_wf%tasks%n_singlet_states = lower_level_n_singlet_states
!
!     Set convergence threshold for lower lying method
!
      ccs_wf%settings%energy_threshold = 1.0D-04 
      ccs_wf%settings%equation_threshold = 1.0D-04 
!
!     Initialize lower level method
!  
      call ccs_wf%init
!
!     Call driver of lower level method
!
      call ccs_wf%drv
!
!     Done with lower lying method
!
      deallocate(ccs_wf)
!
!
!     ::::::::::::::::::::::::::::::::::::::::::::::
!     -::- Construct CNTO transformation matrix -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::
!
!     Open file of CCS solution vectors
!
      call generate_unit_identifier(unit_solution)
!
      open(unit=unit_solution, file='right_eigenvectors', action='read', status='unknown', &
        access='direct', form='unformatted', recl=dp*((wf%n_o)*(wf%n_v)), iostat=ioerror)  
!
      if (ioerror .ne. 0) write(unit_output,*) 'Error while opening solution file', ioerror
!
!     Allocations
!
      call allocator(R_a_i, (wf%n_o)*(wf%n_v), 1)
!
      call allocator(M_i_j, wf%n_o, wf%n_o)
      M_i_j = zero
!
      call allocator(N_a_b, wf%n_v, wf%n_v)
      N_a_b = zero
!
!     Construct M and N
!
      do state = 1, lower_level_n_singlet_states
         read(unit=unit_solution, rec=state) (R_a_i(i , 1), i = 1, (wf%n_o)*(wf%n_v))
!
         call dgemm('T', 'N', &
                     wf%n_o,  &
                     wf%n_o,  &
                     wf%n_v,  &
                     one,     &
                     R_a_i,   &
                     wf%n_v,  &
                     R_a_i,   &
                     wf%n_v,  &
                     one,     &
                     M_i_j,   &
                     wf%n_o)
!
         call dgemm('N', 'T', &
                     wf%n_v,  &
                     wf%n_v,  &
                     wf%n_o,  &
                     one,     &
                     R_a_i,   &
                     wf%n_v,  &
                     R_a_i,   &
                     wf%n_v,  &
                     one,     &
                     N_a_b,   &
                     wf%n_v)
!
      enddo
!
      call deallocator(R_a_i,(wf%n_o)*(wf%n_v), 1)
!
!     Done with file
!
      close(unit_solution, status='delete')
!
!     :: Normalize Tr(N) and Tr(M) ::
!
      trace = 0
!
      do i = 1, wf%n_o
         trace = trace + M_i_j(i, i)
      enddo
!
      call dscal((wf%n_o)*(wf%n_o), one/trace, M_i_j, 1)
!
      trace = 0
!
      do i = 1, wf%n_v
         trace = trace + N_a_b(i, i)
      enddo
!
      call dscal((wf%n_v)*(wf%n_v), one/trace, N_a_b, 1)
!
!     :: Diagonalize M and N matrix ::
!
      call allocator(eigenvalues_o, wf%n_o, 1)
      call allocator(work, 4*(wf%n_o), 1)
      work = zero
!
      call dsyev('V','U', &
                  wf%n_o, &
                  M_i_j, &
                  wf%n_o, &
                  eigenvalues_o, &
                  work, & 
                  4*(wf%n_o), &
                  info)
!
      call deallocator(work, 4*(wf%n_o), 1)
!
      call allocator(eigenvalues_v, wf%n_v, 1)
      call allocator(work, 4*(wf%n_v), 1)
      work = zero
!
      call dsyev('V','U', &
                  wf%n_v, &
                  N_a_b, &
                  wf%n_v, &
                  eigenvalues_v, &
                  work, & 
                  4*(wf%n_v), &
                  info)
!
      call deallocator(work, 4*(wf%n_v), 1)
!
!     :: Reorder M and N ::
!
!     dsyev orderes eigenvalues and corresponding eigenvectors in ascending order.
!     We wish to select active space according to highest eigenvalues, thus we must reorder
!
      call allocator(M, wf%n_o, wf%n_o)
!
      do i = 1, wf%n_o
!
         j = i -1
         M(:,i) = M_i_j(:, wf%n_o - j)
!
      enddo
!
      call deallocator(M_i_j, wf%n_o, wf%n_o)
!
      call allocator(N, wf%n_v, wf%n_v)
!
      do a = 1, wf%n_v
!
         b = a -1
         N(:,a) = N_a_b(:,wf%n_v-b)
!
      enddo
!
      call deallocator(N_a_b, wf%n_v, wf%n_v)
!
!     Transform C to CNTO
!
      call allocator(C_o, wf%n_ao, wf%n_o)
      call allocator(C_v, wf%n_ao, wf%n_v)
      C_o = zero
      C_v = zero      
!
     do i = 1, wf%n_ao
!
        do j = 1, wf%n_o
           ij = index_two(i, j, wf%n_ao)
           C_o(i, j) = wf%mo_coef(ij, 1) 
        enddo
!
        do j = 1, wf%n_v 
           ij = index_two(i, j + wf%n_o, wf%n_ao)
           C_v(i, j) = wf%mo_coef(ij, 1) 
        enddo
!
     enddo
!
      call allocator(C_o_transformed, wf%n_ao, wf%n_o)
      call dgemm('N', 'N',    &
                  wf%n_ao,    &
                  wf%n_o,     &
                  wf%n_o,     &
                  one,        &
                  C_o,        &
                  wf%n_ao,    &
                  M,          &
                  wf%n_o,     &
                  zero,       &
                  C_o_transformed, &
                  wf%n_ao)
!
      call deallocator(C_o, wf%n_ao, wf%n_o)
!
      call allocator(C_v_transformed, wf%n_ao, wf%n_v)
      call dgemm('N', 'N',                   &
                  wf%n_ao,                   &
                  wf%n_v,                    &
                  wf%n_v,                    &
                  one,                       &
                  C_v,                       &
                  wf%n_ao,                   &
                  N,                         &
                  wf%n_v,                    &
                  zero,                      &
                  C_v_transformed,           &
                  wf%n_ao)
!
      call deallocator(C_v, wf%n_ao, wf%n_v)
! 
      do i = 1, wf%n_ao
!

         do j = 1, wf%n_o
            ij = index_two(i, j, wf%n_ao)
            wf%mo_coef(ij, 1) = C_o_transformed(i, j) 
         enddo
!
         do j = 1, wf%n_v 
            ij = index_two(i, j + wf%n_o, wf%n_ao)
            wf%mo_coef(ij, 1) = C_v_transformed(i, j) 
         enddo
!
      enddo
!
      call deallocator(C_o_transformed, wf%n_ao, wf%n_o)
      call deallocator(C_v_transformed, wf%n_ao, wf%n_v)

!
!     :: Determine number of active orbitals ::
!
      sum_o      = 1
      wf%n_CC2_o = 1
!
      do while ((sum_o .gt. wf%mlcc_settings%delta_o) .and. (wf%n_CC2_o .le. wf%n_o))
!
         sum_o = sum_o - eigenvalues_o(wf%n_o - (wf%n_CC2_o - 1), 1)
         wf%n_CC2_o = wf%n_CC2_o + 1
!
      enddo
!
      sum_v      = 1
      wf%n_CC2_v = 1
!
      do while (sum_v .gt. wf%mlcc_settings%delta_v .and. (wf%n_CC2_v .le. wf%n_v))
!
         sum_v = sum_v - eigenvalues_v(wf%n_v - (wf%n_CC2_v - 1), 1)
         wf%n_CC2_v = wf%n_CC2_v + 1
!
      enddo
!
!     Save information to object
!
      wf%first_CC2_o = 1
      wf%first_CC2_v = 1
!
      wf%n_CCS_o = wf%n_o - wf%n_CC2_o
      wf%n_CCS_v = wf%n_v - wf%n_CC2_v
!
      wf%first_CCS_o = 1 + wf%n_CC2_o
      wf%first_CCS_v = 1 + wf%n_CC2_v
!
      call deallocator(eigenvalues_o, wf%n_o, 1)
      call deallocator(eigenvalues_v, wf%n_v, 1)
!
!     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     -::- Finding orbital energies and new block diagonal C matrix -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     Initialize amplitudes and associated attributes
!
      call wf%initialize_amplitudes
      call wf%initialize_omega
!
!     Set the number of parameters in the wavefunction
!     (that are solved for in the ground and excited state solvers) 
!
      wf%n_parameters = wf%n_t1am
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
!
!     Diagonalize active-active block
!
      call allocator(work, 4*(wf%n_CC2_o), 1)
      call allocator(orbital_energies, (wf%n_CC2_o), 1)
      work = zero
!
      call dsyev('V','U',              &
                  (wf%n_CC2_o),        &
                  wf%fock_ij,          &
                  wf%n_o,              &
                  orbital_energies,    &
                  work,                & 
                  4*(wf%n_CC2_o),      &
                  info)
!
      call deallocator(work, 4*(wf%n_CC2_o), 1)
!
      if (info .ne. 0) then
         write(unit_output,*)'WARNING: Diagonalization of active virtual block not successful. '
               stop
      endif
!
!     Setting orbital energies
!
      do j = 1, wf%n_CC2_o
!
         wf%fock_diagonal(j,1) = orbital_energies(j,1)
!
      enddo
!
     call deallocator(orbital_energies, (wf%n_CC2_o), 1)
!
!     Diagonalize inactive-inactive block 
!
      call allocator(work, 4*wf%n_CCS_o, 1)
      call allocator(orbital_energies, wf%n_CCS_o, 1)
      orbital_energies = zero
      work = zero
!
      call dsyev('V','U',                                   &
                  wf%n_CCS_o,                               &
                  wf%fock_ij(wf%first_CCS_o, wf%first_CCS_o), &
                  wf%n_o,                                   &
                  orbital_energies,                         &
                  work,                                     & 
                  4*wf%n_CCS_o,                             &
                  info)
!
      call deallocator(work, 4*wf%n_CCS_o, 1)
!
      if (info .ne. 0) then
         write(unit_output,*)'WARNING: Diagonalization of inactive virtual block not successful.'
               stop
      endif
!
!     Setting orbital energies
!
      do j = 1, wf%n_CCS_o
!
         wf%fock_diagonal(j + wf%first_CCS_o - 1,1) = orbital_energies(j,1)
!
      enddo
      call deallocator(orbital_energies, wf%n_CCS_o, 1)
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
      call allocator(C_o_transformed, wf%n_ao, wf%n_CC2_o)
      call dgemm('N', 'N',    &
                  wf%n_ao,    &
                  wf%n_CC2_o, &
                  wf%n_CC2_o, &
                  one,        &
                  C_o,        &
                  wf%n_ao,    &
                  wf%fock_ij, &
                  wf%n_o,     &
                  zero,       &
                  C_o_transformed, &
                  wf%n_ao)
!
      do i = 1, wf%n_ao
!
         do j = 1, wf%n_CC2_o
            ij = index_two(i, j, wf%n_ao)
            wf%mo_coef(ij, 1) = C_o_transformed(i, j) 
         enddo
!
      enddo
!
      call deallocator(C_o_transformed, wf%n_ao, wf%n_CC2_o)
!
!     Transform C-matrix (inactive occupied block)
!
      call allocator(C_o_transformed, wf%n_ao, wf%n_CCS_o)
      call dgemm('N', 'N',    &
                  wf%n_ao,    &
                  wf%n_CCS_o, &
                  wf%n_CCS_o, &
                  one,        &
                  C_o(1, wf%first_CCS_o),        &
                  wf%n_ao,    &
                  wf%fock_ij(wf%first_CCS_o, wf%first_CCS_o), &
                  wf%n_o,     &
                  zero,       &
                  C_o_transformed, &
                  wf%n_ao)
!
      call deallocator(C_o, wf%n_ao, wf%n_o)
!
      do i = 1, wf%n_ao
!
         do j = 1, wf%n_CCS_o
            ij = index_two(i, j + wf%first_CCS_o -1, wf%n_ao)
            wf%mo_coef(ij, 1) = C_o_transformed(i, j) 
         enddo
!
      enddo
!
      call deallocator(C_o_transformed, wf%n_ao, wf%n_CCS_o)
!
!     :: Vacant orbitals ::
!
!     Diagonalize active-active block
!
      call allocator(work, 4*(wf%n_CC2_v), 1)
      call allocator(orbital_energies, (wf%n_CC2_v), 1)
      work = zero
!
      call dsyev('V','U',              &
                  (wf%n_CC2_v),   &
                  wf%fock_ab,          &
                  wf%n_v,              &
                  orbital_energies,    &
                  work,                & 
                  4*(wf%n_CC2_v), &
                  info)
!
      call deallocator(work, 4*(wf%n_CC2_v), 1)
!
      if (info .ne. 0) then
         write(unit_output,*)'WARNING: Diagonalization of active virtual block not successful. ', orbital_energies
               stop
      endif
!
!     Setting orbital energies
!
      do j = 1, wf%n_CC2_v
!
         wf%fock_diagonal(j + wf%n_o ,1) = orbital_energies(j,1)
!
      enddo
!
     call deallocator(orbital_energies, (wf%n_CC2_v), 1)
!
!     Diagonalize inactive-inactive block 
!
      call allocator(work, 4*wf%n_CCS_v, 1)
      call allocator(orbital_energies, wf%n_CCS_v, 1)
      orbital_energies = zero
      work = zero
!
      call dsyev('V','U',                                   &
                  wf%n_CCS_v,                               &
                  wf%fock_ab(wf%first_CCS_v, wf%first_CCS_v),   &
                  wf%n_v,                                   &
                  orbital_energies,                         &
                  work,                                     & 
                  4*wf%n_CCS_v,                             &
                  info)
!
      call deallocator(work, 4*wf%n_CCS_v, 1)
!
      if (info .ne. 0) then
         write(unit_output,*)'WARNING: Diagonalization of inactive virtual block not successful.'
               stop
      endif
!
!     Setting orbital energies
!
      do j = 1, wf%n_CCS_v
!
         wf%fock_diagonal(j + wf%n_o + wf%first_CCS_v - 1,1) = orbital_energies(j,1)
!
      enddo
      call deallocator(orbital_energies, wf%n_CCS_v, 1)
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
      call allocator(C_v_transformed, wf%n_ao, wf%n_CC2_v)
      call dgemm('N', 'N',    &
                  wf%n_ao,    &
                  wf%n_CC2_v,     &
                  wf%n_CC2_v,     &
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
!
         do j = 1, wf%n_CC2_v
            ij = index_two(i, j + wf%n_o, wf%n_ao)
            wf%mo_coef(ij, 1) = C_v_transformed(i, j) 
         enddo
!
      enddo
!
      call deallocator(C_v_transformed, wf%n_ao, wf%n_CC2_v)
!
!     Transform C-matrix (inactive virtual block)
!
      call allocator(C_v_transformed, wf%n_ao, wf%n_CCS_v)
      call dgemm('N', 'N',    &
                  wf%n_ao,    &
                  wf%n_CCS_v, &
                  wf%n_CCS_v, &
                  one,        &
                  C_v(1, wf%first_CCS_v),        &
                  wf%n_ao,    &
                  wf%fock_ab(wf%first_CCS_v, wf%first_CCS_v), &
                  wf%n_v,     &
                  zero,       &
                  C_v_transformed, &
                  wf%n_ao)
!
      do i = 1, wf%n_ao
!
         do j = 1, wf%n_CCS_v
            ij = index_two(i, j + wf%n_o + wf%first_CCS_v - 1, wf%n_ao)
            wf%mo_coef(ij, 1) = C_v_transformed(i, j) 
         enddo
!
      enddo
!
      call deallocator(C_v_transformed, wf%n_ao, wf%n_CCS_v)
      call deallocator(C_v, wf%n_ao, wf%n_v)
!
   end subroutine cc2_cnto_mlcc2
!
!
   module subroutine print_cnto_info_mlcc2(wf)
!!
!!
      implicit none 
!
      class(mlcc2) :: wf
!
      write(unit_output,'(t3,a40, i3)') 'Number of CC2 occupied orbitals:       ', wf%n_CC2_o
      write(unit_output,'(t3,a40, i3/)')'Number of CC2 virtual orbitals:        ', wf%n_CC2_v
!
      write(unit_output,'(t3,a40, i3)') 'Number of inactive occupied orbitals:  ', wf%n_CCS_o
      write(unit_output,'(t3,a40, i3/)')'Number of inactive virtual orbitals:   ', wf%n_CCS_v
      flush(unit_output)
!
   end subroutine print_cnto_info_mlcc2
end submodule orbital_partitioning