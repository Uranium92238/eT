!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
module sad_tool_class
!
!!
!! SAD tool class
!! Written by Sarai D. Folkestad, Oct 2021
!!
!! Generates SAD densities for a set of atomic centers
!!
!
   use parameters
   use global_out, only: output, timing
!
!
   implicit none
!
   type :: sad_tool
!
      real(dp), private :: threshold
      integer,  private :: max_iterations
!
      character(len=200), private :: ao_density_guess
!
   contains
!
      procedure, public :: generate => generate_sad_tool
!
      procedure, private :: generate_sad_for_center
      procedure, nopass, private :: get_unique_center_indices
!
   end type  sad_tool
!
   interface sad_tool
!
      procedure :: new_sad_tool
!
   end interface sad_tool
!
contains
!
!
   pure function new_sad_tool(threshold, max_iterations) result(this)
!!
!!    New
!!    Written by Sarai D. Folkestad, Oct 2021
!!
      implicit none
!
      real(dp), intent(in) :: threshold
      integer, intent(in), optional :: max_iterations
!
      type(sad_tool) :: this
!
      this%ao_density_guess = 'core'
!
      this%threshold = min(threshold, 1.0d-6)
!
      this%max_iterations = 100
      if (present(max_iterations)) this%max_iterations = max_iterations
!
   end function new_sad_tool
!
   subroutine generate_sad_tool(this, n, centers)
!!
!!    Generate
!!    Written by Tor S. Haugland, Sep 2019
!!
!!    Generates SAD density for every unique (atom, basis) pair.
!!
!!       1: a molecular system with only one atom for every (atom, basis) is created.
!!       2: an UHF wavefunction is created for that system.
!!       3: the wavefunction is run through a SCF solver.
!!       4: the density files are moved to where the wavefunction expects them to be.
!!
!!    Modified by Tor S. Haugland, Nov 2019
!!
!!    Removed deletion of files generated in SCF for SAD to fix a re-occuring bug.
!!    Added "Found SAD". Only loop over unique atoms.
!!
!!    Moved and adapted to SAD generation tool by Sarai D. Folkestad, 2021
!!
      use atomic_center_class,      only: atomic_center
      use timings_class,            only: timings
!
      implicit none
!
      class(sad_tool) :: this
!
      integer :: n
      type(atomic_center), dimension(n) :: centers
!
      integer               :: I
      integer, dimension(n) :: unique_atom_index
!
      type(timings), allocatable :: sad_generation_timer
!
      type(atomic_center)      :: center
!
      sad_generation_timer = timings('SAD generation time', pl='normal')
      call sad_generation_timer%turn_on()
!
      call get_unique_center_indices(n, centers, unique_atom_index)
!
      call timing%mute()
!
      do I = 1, n
!
         if (all(unique_atom_index /= I)) cycle
!
         center = centers(I)
         call this%generate_sad_for_center(center)
!
      enddo
!
      call timing%unmute()
!
      call sad_generation_timer%turn_off()
!
   end subroutine generate_sad_tool
!
!
   subroutine get_unique_center_indices(n, centers, center_indices)
!!
!!    Get unique center indices
!!    Written by Tor S. Haugland and Sarai D. Folkestad, Jun 2021
!!
!!    Collects the indices of unique atoms, excluding ghost atoms.
!!
!
      use string_utilities,    only: index_of_unique_strings
      use atomic_center_class, only: atomic_center
!
      implicit none
!
      integer :: n
      type(atomic_center), dimension(n), intent(in) :: centers
      integer, dimension(n), intent(out) :: center_indices
!
      integer :: n_centers, I
!
      character(len=50), dimension(n) :: atom_and_basis
      integer, dimension(n) :: unique_center_indices
!
      do I = 1, n
!
         atom_and_basis(I) = centers(I)%get_identifier_string()
!
      enddo
!
      call index_of_unique_strings(unique_center_indices, n, atom_and_basis)
!
      center_indices = 0
      n_centers = 0
!
      do I = 1, n
!
         if (all(unique_center_indices /= I)) cycle
!
!        No density to generate for ghost atom
         if (centers(I)%is_ghost()) cycle
!
!        No electrons on this atom
         if (centers(I)%number_ == centers(I)%charge) cycle
!
         n_centers = n_centers + 1
         center_indices(n_centers) = I
!
      enddo
!
   end subroutine get_unique_center_indices
!
!
   subroutine generate_sad_for_center(this, center)
!!
!!    Generate SAD for center
!!    Written by Tor S. Haugland, Sep 2019
!!
      use sequential_file_class,    only: sequential_file
      use atomic_center_class,      only: atomic_center
      use uhf_class,                only: uhf
      use scf_solver_factory_class, only: scf_solver_factory
      use scf_solver_class,         only: scf_solver
!
      implicit none
!
      class(sad_tool),     intent(in) :: this
      type(atomic_center), intent(in) :: center
!
      type(uhf),         allocatable :: sad_wf
      class(scf_solver), allocatable :: sad_solver
!
      type(scf_solver_factory) :: factory
!
      character(len=200)    :: name_
!
      character(len=200)    :: alpha_file_name
      character(len=200)    :: beta_file_name
!
      type(sequential_file) :: alpha_density_file
      type(sequential_file) :: beta_density_file
!
      call output%mute()
!
      sad_wf = uhf(fractional_uniform_valence=.true., &
                   multiplicity=center%get_ground_state_multiplicity())
!
      call sad_wf%prepare([center],  embedding=.false., charge=center%charge)
!
      call sad_wf%set_gradient_threshold(this%threshold)
      call sad_wf%prepare_for_scf(restart=.false., &
                                  skip=.false.,    &
                                  ao_density_guess=this%ao_density_guess)
!
      factory = scf_solver_factory(acceleration_type='none', max_iterations=this%max_iterations)
      call factory%create(sad_wf, sad_solver, restart=.false., skip=.false.)
!
      call sad_solver%run(sad_wf)
!
      call sad_wf%orbital_file%delete_()
      call sad_wf%cleanup()
!
      deallocate(sad_wf)
!
      call output%unmute()
      call output%printf('v', 'Generated atomic density for ' //  &
                         adjustl(center%symbol) // ' using UHF/(a0)', &
                         chars=[center%basis], fs='(t6,a)')
!
!     Move density files to where ao_tool can use them,
!     but first delete SAD if it already exists.
!
      name_ = "sad_" // trim(center%get_identifier_string())
      alpha_file_name = trim(name_) // '_alpha'
      beta_file_name  = trim(name_) // '_beta'
!
      alpha_density_file = sequential_file(alpha_file_name)
      if (alpha_density_file%exists()) call alpha_density_file%delete_()
!
      beta_density_file  = sequential_file(beta_file_name)
      if (beta_density_file%exists())  call beta_density_file%delete_()
!
      alpha_density_file = sequential_file('ao_density_a')
      call alpha_density_file%copy(alpha_file_name)
      call alpha_density_file%delete_()
!
      beta_density_file  = sequential_file('ao_density_b')
      call beta_density_file%copy(beta_file_name)
      call beta_density_file%delete_()
!
   end subroutine generate_sad_for_center
!
!
end module sad_tool_class
