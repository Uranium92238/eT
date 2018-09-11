module atomic_class
!
!!
!!    Atom class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use parameters
   use file_class
   use shell_class
   use basis_set_info
   use memory_manager_class
   use disk_manager_class
!
   implicit none
!
   type :: atomic
!
      character(len=2) :: symbol
      integer(i15)     :: number
!
      character(len=100) :: basis ! name of basis set
!
      integer(i15) :: n_ao ! Number of aos sentered on this atom
!
      integer(i15) :: n_shells
      type(shell), dimension(:), allocatable :: shells
!
      real(dp) :: x
      real(dp) :: y
      real(dp) :: z
!
   contains
!
      procedure          :: set_number       => set_number_atom
      procedure, private :: symbol_to_number => symbol_to_number_atom
!
      procedure, private :: get_n_Aufbau     => get_n_Aufbau_atom
      procedure, private :: get_Aufbau_info  => get_Aufbau_info_atom
!
      procedure :: AD                        => AD_atom
!
      procedure :: read_atomic_density       => read_atomic_density_atomic
!
      procedure :: initialize_shells         => initialize_shells_atomic
      procedure :: destruct_shells           => destruct_shells_atomic
!
      procedure :: cleanup                   => cleanup_atomic
!
   end type atomic
!
contains
!
!
   subroutine set_number_atom(atom)
!!
!!    Set number
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Sets the atomic number from atomic symbol
!!
      implicit none
!
      class(atomic) :: atom
!
      call atom%symbol_to_number()
!
   end subroutine set_number_atom
!
!
   subroutine symbol_to_number_atom(atom)
!!
!!    Symbol to number
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Uses the periodic table to determine atomic number
!!    from atomic symbol
!!
      use periodic_table
!
      implicit none
!
      class(atomic) :: atom
!
      integer(i15) :: i = 0
!
      atom%number = 0
!
      do i = 1, size_periodic_table
!
         if (atomic_symbol(i) == atom%symbol) then
!
            atom%number = i
!
            return
!
         endif
      enddo
!
      if (atom%number == 0) then
!
         write(output%unit,'(/t3,a)') 'Error: illegal atomic symbol, check the eT.inp file '
         stop
!
      endif
!
   end subroutine symbol_to_number_atom
!
!
   subroutine AD_atom(atom, density_diagonal_for_atom)
!!
!!    Atomic density
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Relies on the ordering from libint:
!!       - core and valence (s, p, d, ...) (To be filled)
!!       - polerization (p, d, ...)
!!       - augmentation (s, p, d)
!!
!!    OBS! Note that for orbitals to be filled
!!         all s-orbitals are given first, then all p-orbitals and so on.
!!
!!    To be used for superposition of atomic densities (SOAD)
!!
      implicit none
!
      class(atomic) :: atom
!
      real(dp), dimension(atom%n_ao, 1) :: density_diagonal_for_atom
!
      integer(i15) :: n_Aufbau_shells, Aufbau_to_fill
      integer(i15), dimension(:,:), allocatable :: Aufbau_shell_info
!
      integer(i15) :: i, j, shell, ao_offset, n_electrons
!
!     Find number of sub-shells to fill, allocate array for information on how many
!     instances of a certain sub-shell (depends on the basis set) there are
!     and how many electrons into this sub-shell (depends on angular momentum of sub-shell)
!
      n_Aufbau_shells = atom%get_n_Aufbau()
!
      call mem%alloc_int(Aufbau_shell_info, n_Aufbau_shells, 2) ! [instances of specific Aufbau shell, Order for filling shell]
!
      call atom%get_Aufbau_info(n_Aufbau_shells, Aufbau_shell_info)
!
      shell = 0
!
      density_diagonal_for_atom = zero
      ao_offset = 0
      shell = 0
      n_electrons = atom%number 
!
      do j = 1, n_Aufbau_shells
!
         if (n_electrons == 0) return
!
!        Determine the jth Aufbau shell to fill
!
         do i = 1, n_Aufbau_shells
!
            if (Aufbau_shell_info(i, 2) == j) Aufbau_to_fill = i
!
         enddo
!
!        Determine the shell number just before the first 
!        shell in the Aufbau shell (later on we move do a ++
!        to the first shell)
!
         shell = 0
!
         do i = 1, Aufbau_to_fill - 1
!
            shell = shell + Aufbau_shell_info(i, 1)
!
         enddo 
!
!        Determine the AO offset to just before the first 
!        shell to fill 
!
         ao_offset = 0
!
         do i = 1, shell
!
            ao_offset = ao_offset + atom%shells(i)%size
!
         enddo
!
!        Fill in electron in the correct place 
!
! Eirik stuff - like Sarai, I get junk & it's the same reproduced junk - so: no go
!          shell = shell + 1

! !
!          if (n_electrons .gt. 2*atom%shells(shell)%size) then ! fill entire shell 
! !
!             density_diagonal_for_atom(ao_offset + 1 : ao_offset + atom%shells(shell)%size, 1) &
!                      = two
! !
!          else ! smear out 
! !
!             density_diagonal_for_atom(ao_offset + 1 : ao_offset + atom%shells(shell)%size, 1) &
!                      = real(n_electrons, kind=dp)/real(atom%shells(shell)%size, kind=dp) ! 4/3 4/3 4/3, or n_electrons/size 
! !
!          endif
! End Eirik stuff

         do i = 1, Aufbau_shell_info(Aufbau_to_fill, 1)

            shell = shell + 1

            density_diagonal_for_atom(ao_offset + 1 : ao_offset + atom%shells(shell)%size, 1) &
                = real(min(n_electrons , 2*atom%shells(shell)%size), kind=dp)&
                 /(real(atom%shells(shell)%size, kind=dp))&
                 /(real(Aufbau_shell_info(Aufbau_to_fill, 1), kind=dp))

            ao_offset = ao_offset + atom%shells(shell)%size

         enddo

         n_electrons = n_electrons - min(n_electrons, 2*atom%shells(shell)%size)
!
      enddo
!
      call mem%dealloc_int(Aufbau_shell_info, n_Aufbau_shells, 2)
!
   end subroutine AD_atom
!
!
   integer function get_n_Aufbau_atom(atom)
!!
!!    Get number of Aufbau orbitals
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Returns the number of subshells filled according to Aufbau principle
!!
      implicit none
!
      class(atomic) :: atom
!
      if (atom%number .le. 0) then
!
         write(output%unit, *) 'Error: Illegal atomic number!'
         stop
!
      elseif ((atom%number .gt. 0) .and. (atom%number .le. 2)) then ! H - He
!
         get_n_Aufbau_atom = 1 ! Fill in 1s
!
      elseif ((atom%number .gt. 2) .and. (atom%number .le. 10)) then ! Li - Ne
!
         get_n_Aufbau_atom = 3 ! Fill in 1s, 2s, 2p
!
      elseif ((atom%number .gt. 10) .and. (atom%number .le. 18)) then ! Na - Ar
!
         get_n_Aufbau_atom = 5 ! Fill in 1s, 2s, 3s, 2p, 3p
!
      elseif ((atom%number .gt. 18) .and. (atom%number .le. 20)) then ! K - Ca
!
         get_n_Aufbau_atom = 6 ! Fill in 1s, 2s, 3s, 4s, 2p, 3p!
!
      else
!
         write(output%unit, *) 'Error: eT cannot handle atoms heavier than Ca yet!'
         stop
!
      endif
!
   end function get_n_Aufbau_atom
!
!
   subroutine get_Aufbau_info_atom(atom, n_Aufbau_shells, Aufbau_shell_info)
!!
!!
!!
      implicit none
!
      class(atomic) :: atom
!
      integer(i15) :: n_Aufbau_shells
      integer(i15), dimension(n_Aufbau_shells, 2) :: Aufbau_shell_info
!
      if (atom%number .le. 0) then
!
         write(output%unit, *) 'Error: Illegal atomic number!'
         stop
!
      elseif ((atom%number .gt. 0) .and. (atom%number .le. 2)) then ! H - He
!
         call get_shells_to_fill_H_to_He(atom%basis,  Aufbau_shell_info)
!
      elseif ((atom%number .gt. 2) .and. (atom%number .le. 10)) then ! Li - Ne
!
         call get_shells_to_fill_Li_to_Ne(atom%basis,  Aufbau_shell_info)
!
      elseif ((atom%number .gt. 10) .and. (atom%number .le. 18)) then ! Na - Ar
!
         call get_shells_to_fill_Na_to_Ar(atom%basis,  Aufbau_shell_info)
!
      elseif ((atom%number .gt. 18) .and. (atom%number .le. 20)) then ! K - Ca
!
         call get_shells_to_fill_K_to_Ca(atom%basis,  Aufbau_shell_info)
!
      else
!
         write(output%unit, *) 'Error: eT cannot handle atoms heavier than Ca yet!'
         stop
!
      endif
!
   end subroutine get_Aufbau_info_atom
!
!
   subroutine initialize_shells_atomic(atom)
!!
!!    Initialize shells
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
      implicit none
!
      class(atomic) :: atom
!
      if (.not. allocated(atom%shells)) allocate(atom%shells(atom%n_shells))
!
   end subroutine initialize_shells_atomic
!
!
   subroutine destruct_shells_atomic(atom)
!!
!!    destruct shells
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
      implicit none
!
      class(atomic) :: atom
!
      if (allocated(atom%shells)) deallocate(atom%shells)
!
   end subroutine destruct_shells_atomic
!
!
   subroutine cleanup_atomic(atom)
!!
!!    Cleanup
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
      implicit none
!
      class(atomic) :: atom
!
      call atom%destruct_shells
!
   end subroutine cleanup_atomic
!
!
   subroutine read_atomic_density_atomic(atom, atomic_density)
!!
!!    Read atomic density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Read the atomic density matrices from file and adds them 
!!    together. By assumption, these densities are the result of 
!!    atomic UHF calculations with the restriction that valence 
!!    electrons are smeared out to ensure that the density is 
!!    spherically symmetric (as required for a rotationally
!!    invariant SAD guess). 
!!
      implicit none 
!  
      class(atomic) :: atom 
!
      real(dp), dimension(atom%n_ao, atom%n_ao) :: atomic_density
!
      real(dp), dimension(:,:), allocatable :: temporary
      character(len=100)                    :: alpha_fname, beta_fname
!
      character(len=255) :: sad_directory
!
      type(file) :: alpha_density
      type(file) :: beta_density
!
      integer(i15) :: i, j, ioerror
!
      atomic_density = zero   
!
      call get_environment_variable("SAD_ET_DIR", sad_directory)
      write(output%unit, *)sad_directory
!
      alpha_fname = trim(sad_directory) // '/' // trim(atom%basis) // '/' // trim(atom%symbol) // '_' // 'alpha' // '.inp'
      beta_fname  = trim(sad_directory) // '/' //trim(atom%basis) // '/' // trim(atom%symbol) // '_' // 'beta' // '.inp'

   !  alpha_fname = '/Users/sarai/new_eT/src/molecular_system/sad' &
   !            // '/' // trim(atom%basis) // '/' // trim(atom%symbol) // '_' // 'alpha' // '.inp'
   !  beta_fname  = '/Users/sarai/new_eT/src/molecular_system/sad' &
   !            // '/' //trim(atom%basis) // '/' // trim(atom%symbol) // '_' // 'beta' // '.inp'

!
      call alpha_density%init(trim(alpha_fname), 'sequential', 'formatted')
      call beta_density%init(trim(beta_fname), 'sequential', 'formatted')
!
      call disk%open_file(alpha_density, 'read', 'rewind')
      call disk%open_file(beta_density, 'read', 'rewind')
!
      call mem%alloc(temporary, atom%n_ao, atom%n_ao)
      read(alpha_density%unit, *, iostat=ioerror) ((temporary(i,j), j = 1, atom%n_ao), i = 1, atom%n_ao)
      atomic_density = atomic_density + temporary
!
      read(beta_density%unit, *, iostat=ioerror) ((temporary(i,j), j = 1, atom%n_ao), i = 1, atom%n_ao)
      atomic_density = atomic_density + temporary
!
      call mem%dealloc(temporary, atom%n_ao, atom%n_ao)
!
      call disk%close_file(alpha_density)
      call disk%close_file(beta_density)
!
   end subroutine read_atomic_density_atomic
!
!
end module atomic_class
