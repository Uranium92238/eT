submodule (mlcc2_class) orbital_partitioning
!
!!
!!    Orbital partitioning submodule (MLCC2) 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2017
!!
!!    Contains the following family of procedures of the MLCC2 class:
!!
!!    orbital_partitioning:   Directs the orbital partitioning
!!    cholesky_localization:  Directs orbital localization by cholesky decomposition  
!!    cholesky_orbital_drv:   Directs construction of new orbitals                     
!!    cholesky_decomposition: Cholesky decomposes the density
!!    cholesky_orbitals:      Constructs new orbitals (C matrix) from cholesky vectors 
!! 
!!    Contains the following subroutines and functions
!
   implicit none 
!
   logical :: debug = .true.
!
!
contains
!
!
   module subroutine orbital_partitioning_mlcc2(wf)
!!
!!
      implicit none
!
      class(mlcc2) :: wf
!
      if (wf%mlcc_settings%cholesky) then

         call wf%cholesky_localization
!
      endif
!
   end subroutine orbital_partitioning_mlcc2
!
!
   module subroutine cholesky_localization_mlcc2(wf)
!!
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(:,:), allocatable     :: orbitals
      real(dp), dimension(:,:), allocatable     :: orbital_energies
      integer(i15)                              :: n_nuclei
      integer(i15), dimension(:,:), allocatable :: ao_center_info, n_ao_on_center
      integer(i15)                              :: active_space = 0, n_CC2_atoms = 0
      integer(i15)                              :: offset_o = 0, offset_v = 0
      integer(i15), dimension(:,:), allocatable :: active_atoms
      integer(i15), dimension(:,:), allocatable :: active_ao_index_list
      integer(i15)                              :: n_active_aos = 0
      integer(i15), dimension(:,:), allocatable :: n_vectors_o, n_vectors_v
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
      real(dp), dimension(:,:), allocatable :: C
!
!     Start timings
!
      call cpu_time(start_chol_deco)
!
!     :: Prepare for decomposition ::
!
!     Prints
!
      write(unit_output,'(/t3,a/)')    ':: Cholesky decomposition '
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
!     Cholesky localized orbitals
!
      call allocator (orbitals, wf%n_ao, wf%n_mo)
      call allocator (orbital_energies, wf%n_mo, 1)
!
!
!     Print number of active spaces requested active spaces
!
      write(unit_output,'(t3,a,i3,a/)')'Requested ', wf%n_active_spaces,' active space(s).'
!
!     Construct canonical occupied and vacant density matrices     
!
      call allocator(density_o, wf%n_ao, wf%n_ao)
      call allocator(density_v, wf%n_ao, wf%n_ao)
      density_o = zero
      density_v = zero
!
      call wf%construct_density_matrices(density_o, density_v)
!
!     Construct ao fock matrix
!
      call allocator(ao_fock, wf%n_ao, wf%n_ao)
      call wf%construct_ao_fock(ao_fock, density_o)
!
!     Prepare for loop over active spaces
!
      write(unit_output,'(t3,a,5x,a,5x,a)')'Active space:','Number of atoms:','Atoms:'
!
      offset_o = 1
      offset_v = 1 + wf%n_o
!
!     Variables for storing information on spaces
!
      call allocator_int(n_vectors_o, wf%n_active_spaces + 1, 1)
      call allocator_int(n_vectors_v, wf%n_active_spaces + 1, 1)
!
!     Start loop over active spaces
!
      do active_space = 1, wf%n_active_spaces
!
!        Get CC2-active atoms
!
         n_CC2_atoms =  get_number_of_active_atoms(unit_cholesky_decomp, active_space, 'CC2  ')
!
!        Allocate active atoms list
!
         call allocator_int(active_atoms, n_CC2_atoms, 1)
!
!        Get list of active atoms
!
         call get_active_atoms(unit_cholesky_decomp, active_atoms,  n_CC2_atoms, active_space, 'CC2  ')
!
!        Sanity check on active atoms
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
!        Prints
!
         write(unit_output,'(t3,i3,15x,i3,18x,i3,i3)')active_space, n_CC2_atoms, active_atoms
         flush(unit_output)
!
!        :: Decomposing the occupied density ::
!
         n_active_aos = 0
         do i = 1, n_CC2_atoms
            n_active_aos = n_active_aos + n_ao_on_center(active_atoms(i,1),1)
         enddo
!
!        Construct active_ao_index_list
!
         call allocator_int(active_ao_index_list, n_active_aos, 1)
         call construct_active_ao_index_list(active_ao_index_list, n_active_aos, active_atoms, &
                                                  n_CC2_atoms, ao_center_info, wf%n_ao)
!
!        :: Occupied part ::
!
         call wf%cholesky_orbital_drv(orbitals, orbital_energies, offset_o, ao_fock, density_o, n_vectors_o(active_space,1),&
                              .true., n_active_aos, active_ao_index_list)
!
!        :: Virtual part ::
!
         call wf%cholesky_orbital_drv(orbitals, orbital_energies, offset_v, ao_fock, density_v, n_vectors_v(active_space,1),&
                              .true., n_active_aos, active_ao_index_list)
!
!        Calculate new offset         
!
         offset_o = offset_o + n_vectors_o(active_space, 1)
         offset_v = offset_v + n_vectors_v(active_space, 1)
!
         call deallocator_int(active_atoms, n_CC2_atoms, 1)
         call deallocator_int(active_ao_index_list, n_active_aos, 1)
!
!        Save active space information
!
         wf%n_CC2_o(active_space, 1) = n_vectors_o(active_space, 1)
         wf%n_CC2_v(active_space, 1) = n_vectors_v(active_space, 1)
!   
      enddo
!
!     :: CCS orbitals ::
!
!     Occupied part ::
!
      call wf%cholesky_orbital_drv(orbitals, orbital_energies, offset_o, ao_fock, density_o, n_vectors_o(wf%n_active_spaces+1,1),&
                              .false., n_active_aos)
!
!     :: Virtual part ::
!
      call wf%cholesky_orbital_drv(orbitals, orbital_energies, offset_v, ao_fock, density_v, n_vectors_v(wf%n_active_spaces+1,1),&
                              .false., n_active_aos)
!
      call deallocator(density_v, wf%n_ao, wf%n_ao)   
      call deallocator(density_o, wf%n_ao, wf%n_ao)   
      call deallocator(ao_fock, wf%n_ao, wf%n_ao)
!
!     Save inactive space information
!
      wf%n_CCS_o = n_vectors_o(wf%n_active_spaces + 1, 1)
      wf%n_CCS_v = n_vectors_v(wf%n_active_spaces + 1, 1)
!
      call deallocator_int(n_vectors_o, wf%n_active_spaces + 1, 1)
      call deallocator_int(n_vectors_v, wf%n_active_spaces + 1, 1)
!
!     Check orthogonality
!
!     NEED OVERLAP MATRIX!
!
      wf%mo_coef       = orbitals
      wf%fock_diagonal = orbital_energies
!
      call deallocator (orbitals, wf%n_ao, wf%n_mo)
      call deallocator (orbital_energies, wf%n_mo, 1)
      call deallocator_int(n_ao_on_center, n_nuclei, 2)
      call deallocator_int(ao_center_info, wf%n_ao, 2)
!     
!     Close cholesky.inp
!
      close(unit_cholesky_decomp)
!
!     Print timings
!
      call cpu_time(end_chol_deco)
      write(unit_output,'(/t3,a27,f14.8/)') 'Total time (seconds):', end_chol_deco - start_chol_deco
      flush(unit_output)
!
!
   end subroutine cholesky_localization_mlcc2
!
!
   module subroutine cholesky_orbital_drv_mlcc2(wf, orbitals, orbital_energies, offset, ao_fock, density, n_vectors,&
                              selection, n_active_aos, active_ao_index_list)
!!
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
!     Construct cholesky vectors
!
      if (selection) then
         call wf%cholesky_decomposition(density, cholesky, &
                                               n_vectors, .true., n_active_aos, active_ao_index_list)
      else 
         call wf%cholesky_decomposition(density, cholesky, &
                                               n_vectors, .false., n_active_aos)
      endif   
!
!     Construct cholesky orbitals
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
   end subroutine cholesky_orbital_drv_mlcc2
!
!
   module subroutine cholesky_decomposition_mlcc2(wf, density, cholesky_vectors,&
                                                     n_vectors, selection, n_active_aos, active_ao_index_list)
!!
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
!!
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
!     :: Construct ao-fock ::
!
!
     if (n_vectors == 0) return
!
      call allocator(X, n_vectors, wf%n_ao)
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

!
      call deallocator(X, n_vectors, wf%n_ao)
!
!     Diagonalize mo-fock 
!
      call allocator(work, 4*n_vectors, 1)
      work = zero
!
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
!
!
      call deallocator(work, 4*n_vectors, 1)
!
      if (info .ne. 0) then
         write(unit_output,*)'WARNING: Diagonalization not successful.'
               stop
      endif
!
!     Construct orbitals
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
   module function get_number_of_active_spaces(unit_cholesky_decomp)
!!
!!    Get number of active spaces
!!    Written by Sarai D. Folkestad June 2017
!!
!!    Reads cholesky.inp, and returns the number of active spaces for the mlcc calculation.
!!
      implicit none
!
      integer(i15) :: get_number_of_active_spaces
      integer(i15) :: unit_cholesky_decomp
!
      character(len=40) :: line
      integer(i15)      :: ioerror = 0
!
      do
         read(unit_cholesky_decomp,'(a40)')line
!
         do while (line(1:1) == '!' .or. trim(line) == '') ! Comment or blank line: read the next line
!
            read(unit_cholesky_decomp,'(a40)') line
         enddo
!
         if (trim(line) == 'Active spaces:') then
!
            read(unit_cholesky_decomp,'(i5)', iostat=ioerror) get_number_of_active_spaces
            if (ioerror .ne. 0) then
               write(unit_output,*)'WARNING: Error while reading number of active spaces from cholesky.inp'
               stop
            endif
            exit
!
         elseif(trim(line) == '#end of Cholesky input') then
!
            backspace(unit_cholesky_decomp)
            write(unit_output,*)'WARNING: MLCC wavefunction requested, but number of active spaces not given in cholesky.inp'
            stop
!
         endif
      enddo
!
   end function get_number_of_active_spaces
!
   module function get_number_of_active_atoms(unit_cholesky_decomp, active_space, ml_level)
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
      integer(i15)      :: active_space         ! Which active space is considered 
      character(len=5)  :: ml_level             ! CC2/CCSD/CC3
!
      character(len=40) :: line
      integer(i15)      :: ioerror = 0
      integer(i15)      :: previous_spaces, previous_atoms, previous_n_atoms
!
      rewind(unit_cholesky_decomp)
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
            do previous_spaces = 1, active_space -1
               read(unit_cholesky_decomp,'(i5)', iostat=ioerror) previous_n_atoms
               read(unit_cholesky_decomp,*, iostat=ioerror) previous_atoms
            enddo
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
   module subroutine get_active_atoms(unit_cholesky_decomp, active_atoms, n_active_atoms, active_space, ml_level)
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
      integer(i15)      :: active_space, n_active_atoms
      character(len=5)  :: ml_level
!
      integer(i15), dimension(n_active_atoms,1) :: active_atoms
!
      character(len=40) :: line
      integer(i15)      :: ioerror = 0
      integer(i15)      :: previous_spaces, previous_atoms, previous_n_atoms
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
            do previous_spaces = 1, active_space -1
               read(unit_cholesky_decomp,*, iostat=ioerror) previous_n_atoms
               read(unit_cholesky_decomp,*, iostat=ioerror) previous_atoms
            enddo
            read(unit_cholesky_decomp,*, iostat=ioerror) previous_n_atoms
!
            read(unit_cholesky_decomp,*, iostat=ioerror) active_atoms
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
!!
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
!!
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
!
end submodule orbital_partitioning