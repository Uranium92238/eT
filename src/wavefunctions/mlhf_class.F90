module mlhf_class
!
!!
!!    Multilevel Hartree-Fock (MLHF) class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use hf_class
!
   use reordering
   use interval_class
   use index
   use active_atoms_class
!
   use array_utilities
!
   implicit none
!
   type, extends(hf):: mlhf
!
      type(active_atoms) :: active_space
!
   contains
!
      procedure :: initialize => initialize_mlhf
      procedure :: finalize => finalize_mlhf
!
      procedure :: read_info => read_info_mlhf
!
      procedure :: construct_virtual_density => construct_virtual_density_mlhf
!
   end type mlhf
!
!
contains
!
!
   subroutine initialize_mlhf(wf)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(mlhf) :: wf
!
      wf%name = 'MLHF'
!
      call wf%read_info
!
   end subroutine initialize_mlhf
!
!
   subroutine finalize_mlhf(wf)
!!
!!    Finalize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(mlhf) :: wf
!
!     Nothing here yet
!
   end subroutine finalize_mlhf
!
!
   subroutine decompose_density_active_mlhf(wf)
!
      implicit none
!
      class(mlhf) :: wf
!
      real(dp), dimension(:,:), allocatable :: cholesky_vectors
!
      integer(i15):: i, j, k, n_active_aos, ao_offset, active_ao_counter, n_vectors
!
      integer(i15), dimension(:,:), allocatable :: active_aos

!
      n_active_aos = 0
!
      do i = 1, wf%active_space%n_active_atoms
!
         n_active_aos = n_active_aos + wf%system%atoms(wf%active_space%atoms(i, 1))%n_ao
!
      enddo 
!
      call mem%alloc_int(active_aos, n_active_aos, 1)
!
      ao_offset = 0
      active_ao_counter = 0
!
      do i = 1, wf%system%n_atoms
         do j = 1, wf%active_space%n_active_atoms
!
            if (i == wf%active_space%atoms(j, 1)) then
!
               do k = 1, wf%system%atoms(i)%n_ao
!
                  active_ao_counter = active_ao_counter + 1
!
                  active_aos(active_ao_counter, 1) = k + ao_offset
!   
               enddo
!
            endif
!
         enddo
!
         ao_offset = ao_offset + wf%system%atoms(i)%n_ao
!
      enddo 
!
      call mem%alloc(cholesky_vectors, wf%n_ao, n_active_aos)
!
      call cholesky_decomposition_limited_diagonal(wf%ao_density, cholesky_vectors, wf%n_ao, &
                                                     n_vectors, 1.0d-9, n_active_aos, active_aos)
!
      call mem%dealloc_int(active_aos, n_active_aos, 1)
      call mem%dealloc(cholesky_vectors, wf%n_ao, n_active_aos)
!
   end subroutine decompose_density_active_mlhf
!
!
   subroutine read_info_mlhf(wf)
!
!!
!!    Read information
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!     mlhf
!!        n_active_atoms
!!        active_atoms
!!     end mlhf
!!
      implicit none
!
      class(mlhf) :: wf
!
      character(len=100) :: line
      character(len=100) :: current_basis
!
      integer(i15) :: i = 0,j = 0, ioerror, first, last
!
      call input%init('eT.inp', 'sequential', 'formatted')
      call disk%open_file(input, 'read')
      rewind(input%unit)
!
      read(input%unit,'(a)', iostat=ioerror) line
      line = remove_preceding_blanks(line)
!
      do while ((trim(line) .ne. 'end mlhf') .and. (ioerror .eq. 0))
!
         if (trim(line) == 'mlhf') then ! found cholesky section in input
!
            read(input%unit, *) wf%active_space%n_active_atoms
!
            call mem%alloc_int(wf%active_space%atoms, wf%active_space%n_active_atoms, 1)
!
            read(input%unit,'(a)', iostat=ioerror) line
            line = remove_preceding_blanks(line)
!
            if (line(1:1)=='[') then ! range given
!
               do i = 2, 100
!
                  if (line(i:i) == ',') exit
!
               enddo
!
               read(line(2:i-1), *) first
!
               do j = i, 100
!
                  if (line(j:j) == ']') exit
!
               enddo
!
               read(line(i+1:j-1), *) last
!
               if ((last - first + 1) .ne. wf%active_space%n_active_atoms) then
!
                  write(output%unit, *) 'Error: actual number of active atoms does not match number provided in input.'
                  stop
!
               endif
!
               do i = first, last
!
                  wf%active_space%atoms(i - first + 1, 1) = i
!
               enddo
!
            else ! specific atoms are given
!
               read(input%unit, *) wf%active_space%atoms
!
            endif
!
         endif
!
         read(input%unit,'(a)') line
         line = remove_preceding_blanks(line)
!
      enddo
      backspace(input%unit)
!
      call disk%close_file(input)
!
   end subroutine read_info_mlhf
!
!
  subroutine construct_virtual_density_mlhf(wf)
!!
!!
      implicit none 
!
      class(mlhf) :: wf
!
  end subroutine construct_virtual_density_mlhf
!
!  1. Hvordan skal jeg få rutiner fra HF (midlertidig arve derfra?)
!  2. Skrive inputleser
!  3. Teste dekomponering
!  4. Lage V_xy
!  5. ML
!
!
end module mlhf_class
