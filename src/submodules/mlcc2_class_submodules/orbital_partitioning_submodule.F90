submodule (mlcc2_class) orbital_partitioning
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
!
      call wf%cholesky_orbital_localization
!
      endif
!
   end subroutine orbital_partitioning_mlcc2
!
!
   module subroutine cholesky_orbital_localization_mlcc2(wf)
!!
!!
      implicit none
!
      class(mlcc2) :: wf
!
      integer(i15) :: n_nuclei = 0, nucleus = 0, active_space = 0, n_CC2_atoms = 0
!
      integer(i15) :: unit_center = 0, unit_cholesky_decomp = 0
      integer(i15) :: ioerror     = 0
      integer(i15) :: offset      = 0
      integer(i15) :: i           = 0
!
      logical :: file_exists
!
      integer(i15), dimension(:,:), allocatable :: ao_center_info, n_ao_on_center
      integer(i15), dimension(:,:), allocatable :: active_atoms
      integer(i15), dimension(:,:), allocatable :: active_ao_index_list
!
      real(dp), dimension(:,:), allocatable :: density_o
      real(dp), dimension(:,:), allocatable :: density_v
      real(dp), dimension(:,:), allocatable :: cholesky_vectors_o
      real(dp), dimension(:,:), allocatable :: cholesky_vectors_v
!
      integer(i15) :: n_active_aos = 0
      integer(i15) :: n_vectors_o  = 0
      integer(i15) :: n_vectors_v  = 0
!
!     :: Read center_info ::
!
      call generate_unit_identifier(unit_center)
      open(unit=unit_center, file='center_info', status='unknown', form='unformatted', iostat=ioerror)
      if (ioerror .ne. 0) write(unit_output)'WARNING: Error while opening center_info'
      rewind(unit_center)
!
      read(unit_center) n_nuclei, wf%n_ao
!
      call allocator_int(n_ao_on_center, n_nuclei, 1)
!
      read(unit_center, iostat=ioerror) n_ao_on_center
      if (ioerror .ne. 0) write(unit_output,*)'WARNING: Error while reading center_info'
!
      call allocator_int(ao_center_info, wf%n_ao, 2)
!
      offset = 1
      do nucleus = 1, n_nuclei
!
         do i = 1, n_ao_on_center(nucleus, 1)
            ao_center_info(offset + i - 1,1)=nucleus
         enddo
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
!     :: Prepare for decomposition ::
!
!     Prints
!
      write(unit_output,'(/t3,a/)')    ':: Cholesky decomposition '
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
      if (ioerror .ne. 0) write(unit_output)'WARNING: Error while opening cholesky.inp'
      rewind(unit_cholesky_decomp)
!
!     Get number of active spaces
!
      wf%n_active_spaces = get_number_of_active_spaces(unit_cholesky_decomp)
      write(unit_output,'(t3,a,i3,a/)')'Requested ', wf%n_active_spaces,' active space(s).'
!
!     Construct canonical occupied density matrix     
!
      call allocator(density_o, wf%n_ao, wf%n_ao)
      density_o = zero
!
      call dgemm('N', 'T',    &
                  wf%n_ao,    &
                  wf%n_ao,    &
                  wf%n_o,     &
                  one,        &
                  wf%mo_coef, &
                  wf%n_ao,    &
                  wf%mo_coef, &
                  wf%n_ao,    &
                  zero,       &
                  density_o,  &
                  wf%n_ao)
!
      call allocator(density_v, wf%n_ao, wf%n_ao)
      density_v = zero
!
      call dgemm('N', 'T',                   &
                  wf%n_ao,                   &
                  wf%n_ao,                   &
                  wf%n_v,                    &
                  one,                       &
                  wf%mo_coef(1, wf%n_o + 1), &
                  wf%n_ao,                   &
                  wf%mo_coef(1, wf%n_o + 1), &
                  wf%n_ao,                   &
                  zero,                      &
                  density_o,                 &
                  wf%n_ao)
!
!     Loop over active spaces
!
      write(unit_output,'(t3,a,5x,a,5x,a)')'Active space:','Number of atoms:','Atoms:'
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
         do i = 1, n_active_atoms
            if (active_atoms(i,1) .gt. n_nuclei) then
               write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
               stop
            elseif (n_active_atoms .gt. n_nuclei) then
               write(unit_output,*) 'WARNING: Illegal chioce of atoms in cholesky.inp.'
               stop
            endif
         enddo
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
         call allocator(cholesky_vectors_o, wf%n_ao, wf%n_ao)
         cholesky_vectors_o = zero
         call allocator(cholesky_vectors_v, wf%n_ao, wf%n_ao)
         cholesky_vectors_v = zero
!
         call wf%cholesky_decomposition(density_o, cholesky_vectors_o, n_active_aos, active_ao_index_list, n_vectors_o)
         call wf%cholesky_decomposition(density_v, cholesky_vectors_v, n_active_aos, active_ao_index_list, n_vectors_v)
!
!        Construct new active orbitals
!
         call deallocator(cholesky_vectors_o, wf%n_ao, wf%n_ao)
         call deallocator(cholesky_vectors_v, wf%n_ao, wf%n_ao)
!
         call deallocator_int(active_atoms, n_CC2_atoms, 1)
         call deallocator_int(active_ao_index_list, n_active_aos, 1)

      enddo
!
      call deallocator(density_o, wf%n_ao, wf%n_ao) 
      call deallocator(density_v, wf%n_ao, wf%n_ao) 
!
!     Close cholesky.inp
!
      close(unit_cholesky_decomp)
      call deallocator_int(n_ao_on_center, n_nuclei, 2)
      call deallocator_int(ao_center_info, wf%n_ao, 2)
!
   end subroutine cholesky_orbital_localization_mlcc2
!
!
   module subroutine cholesky_decomposition_mlcc2(wf, density, cholesky_vectors, n_active_aos, active_ao_index_list, n_vectors)
!!
!!
      implicit none
!
      class(mlcc2)                                   :: wf
      integer(i15)                                   :: n_active_aos
      integer(i15)                                   :: n_vectors
      integer(i15), dimension( n_active_aos,1)       :: active_ao_index_list
      real(dp), dimension(wf%n_ao,wf%n_ao)           :: density
      real(dp), dimension(wf%n_ao, wf%n_ao)          :: cholesky_vectors
!
      integer(i15) :: index_max = 0
      real(dp)     :: max_diagonal
!
      real(dp), parameter  :: threshold = 1.0d-07
      real(dp), parameter  :: tolerance = 1.0d-09
!
      integer(i15) :: i = 0, j = 0, k = 0
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
   module subroutine cholesky_orbitals_mlcc2(wf, cholesky_vectors, n_vectors, orbitals)
!!
!!
!!
      implicit none
!
      class(mlcc2)                            :: wf
      real(dp), dimension(wf%n_ao, wf%n_ao)   :: cholesky_vectors
      real(dp), dimension(wf%n_ao, n_vectors) :: orbitals
      integer(i15)                            :: n_vectors
!
!     Construct ao-fock  
!
!
!     Transform ao-fock
!
!
!     Diagonalize ao-fock 
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
end submodule orbital_partitioning