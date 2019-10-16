!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
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
module atomic_class
!
!!
!!    Atom class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use parameters
   use global_out, only : output
   use sequential_file_class, only : sequential_file
   use shell_class, only : shell
   use memory_manager_class, only : mem
!
   implicit none
!
   type :: atomic
!
      character(len=2) :: symbol
      integer     :: number_
!
      character(len=100) :: basis ! name of basis set
!
      integer :: n_ao ! Number of aos sentered on this atom
!
      integer :: n_shells
      type(shell), dimension(:), allocatable :: shells
!
      real(dp) :: x
      real(dp) :: y
      real(dp) :: z
!
      integer :: input_number ! Index of atom on input
!
   contains
!
      procedure          :: set_number       => set_number_atom
      procedure, private :: symbol_to_number => symbol_to_number_atom
!
      procedure :: read_atomic_density       => read_atomic_density_atomic
!
      procedure :: initialize_shells         => initialize_shells_atomic
      procedure :: destruct_shells           => destruct_shells_atomic
!
      procedure :: cleanup                   => cleanup_atomic
!
      procedure :: get_multiplicity     => get_multiplicity_atomic
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
      integer :: i = 0
!
      atom%number_ = 0
!
      do i = 1, size_periodic_table
!
         if (atomic_symbol(i) == atom%symbol) then
!
            atom%number_ = i
!
            return
!
         endif
      enddo
!
      if (atom%number_ == 0) then
!
         call output%error_msg('illegal atomic symbol, check the eT.inp file ')
!
      endif
!
   end subroutine symbol_to_number_atom
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
      type(sequential_file) :: alpha_density_file
      type(sequential_file) :: beta_density_file
!
      atomic_density = zero   
!
      alpha_fname = 'sad_' // trim(atom%basis) // '_' // &
                  & trim(atom%symbol) // '_' // 'alpha'
      beta_fname  = 'sad_' // trim(atom%basis) // '_' // &
                  & trim(atom%symbol) // '_' // 'beta'
!
      alpha_density_file = sequential_file(trim(alpha_fname))
      beta_density_file  = sequential_file(trim(beta_fname))
!
      call alpha_density_file%open_('read', 'rewind')
      call beta_density_file%open_( 'read', 'rewind')
!
      call mem%alloc(temporary, atom%n_ao, atom%n_ao)
!
      call alpha_density_file%read_(temporary, atom%n_ao*atom%n_ao)
      atomic_density = atomic_density + temporary
!
      call beta_density_file%read_(temporary, atom%n_ao*atom%n_ao)
      atomic_density = atomic_density + temporary
!
      call mem%dealloc(temporary, atom%n_ao, atom%n_ao)
!
      call alpha_density_file%close_
      call beta_density_file%close_
!
   end subroutine read_atomic_density_atomic
!
!
   integer function get_multiplicity_atomic(atom) result(multiplicity)
!!
!!    Get multiplicity
!!    Written by Tor S. Haugland, 2019
!!
!!    Returns the multiplicity from the atomic number.
!!    Based on the table 5.1 in Griffith's "Introduction to Quantum Mechanics", 2. ed
!!
      implicit none
      class(atomic) :: atom
      integer, dimension(36) :: multiplicity_list
!
      if (atom%number_ > 36) call output%error_msg("Atom multiplicity not supported for Z > 36")
!
      multiplicity_list = (/2,1,2,1,2,3,4,3,2,1, &
                            2,1,2,3,4,3,2,1,2,1, &
                            2,3,4,7,6,5,4,3,2,1, &
                            2,3,4,3,2,1          &
                              /)
!
      multiplicity = multiplicity_list(atom%number_)
!
   end function get_multiplicity_atomic
!
!
end module atomic_class
