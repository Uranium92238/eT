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
module hf_visualization_task_class
!
!!
!! HF visualization task class
!! Written by Alexander C. Paul, May 2022
!!
!
   use parameters
   use hf_class,             only: hf
   use hf_task_class,        only: hf_task
   use visualization_class,  only: visualization
   use memory_manager_class, only: mem
!
   implicit none
!
   type, extends(hf_task) :: hf_visualization_task
!
      logical, private  :: plot_orbitals, plot_density, plot_active_density
!
      class(visualization), allocatable, private :: visualizer
!
   contains
!
      procedure, public :: execute &
                        => execute_hf_visualization_task
!
      procedure, private :: visualize_orbitals
      procedure, private :: visualize_density
      procedure, private :: visualize_active_density
!
   end type hf_visualization_task
!
!
   interface hf_visualization_task
!
      procedure :: new_hf_visualization_task
!
   end interface hf_visualization_task
!
!
contains
!
!
   function new_hf_visualization_task() result(this)
!!
!!    New HF visualization task
!!    Written by Alexander C. Paul, May 2022
!!
      use global_in, only: input
!
      implicit none
!
      type(hf_visualization_task) :: this
!
      this%name_ = 'Plotting orbitals and/or the HF density'
!
      this%plot_orbitals = input%is_keyword_present('plot hf orbitals', 'visualization')
      this%plot_density = input%is_keyword_present('plot hf density', 'visualization')
      this%plot_active_density = input%is_keyword_present('plot hf active density', &
                                                          'visualization')
!
   end function new_hf_visualization_task
!
!
   subroutine execute_hf_visualization_task(this, wf)
!!
!!    Execute
!!    Written by Alexander C. Paul, May 2022
!!
      implicit none
!
      class(hf_visualization_task), intent(inout) :: this
!
      class(hf), intent(inout), target :: wf
!
      logical :: plot_anything
!
      plot_anything = any([this%plot_orbitals, this%plot_density, this%plot_active_density])
!
      if (.not. plot_anything) return
!
      call this%print_header()
      call this%start_timer()
!
      this%visualizer = visualization(wf%ao)
      call this%visualizer%initialize(wf%ao)
!
      call this%visualize_orbitals(wf)
      call this%visualize_density(wf)
      call this%visualize_active_density(wf)
!
      call this%visualizer%cleanup()
!
      call this%end_timer()
!
   end subroutine execute_hf_visualization_task
!
!
   subroutine visualize_orbitals(this, wf)
!!
!!    Visualize orbitals
!!    Written by Sarai D. Folkestad, Oct 2019
!!
      use global_in, only: input
!
      implicit none
!
      class(hf_visualization_task), intent(inout) :: this
      class(hf), intent(in) :: wf
!
      integer :: i, n_orbitals_to_plot
!
      integer, dimension(:), allocatable :: orbitals_to_plot
!
      real(dp), dimension(:,:), allocatable :: coefficients
!
      character(len=200), dimension(:), allocatable :: orbital_file_tags
!
      if (.not. this%plot_orbitals) return
!
      n_orbitals_to_plot = input%get_n_elements_for_keyword('plot hf orbitals', 'visualization')
!
      call mem%alloc(orbitals_to_plot, n_orbitals_to_plot)
      call input%get_array_for_keyword('plot hf orbitals', 'visualization', &
                                        n_orbitals_to_plot, orbitals_to_plot)
!
      call mem%alloc(coefficients, wf%ao%n, n_orbitals_to_plot)
      allocate(orbital_file_tags(n_orbitals_to_plot))
!
      do i = 1, n_orbitals_to_plot
!
         call dcopy(wf%ao%n, wf%orbital_coefficients(1, orbitals_to_plot(i)), 1, &
                    coefficients(1, i), 1)
!
         write(orbital_file_tags(i), '(a, i4.4)') 'MO_', orbitals_to_plot(i)
!
      enddo
!
      call mem%dealloc(orbitals_to_plot, n_orbitals_to_plot)
!
      call this%visualizer%plot_orbitals(wf%ao, coefficients, &
                                         n_orbitals_to_plot, orbital_file_tags)
!
      call mem%dealloc(coefficients, wf%ao%n, n_orbitals_to_plot)
      deallocate(orbital_file_tags)
!
   end subroutine visualize_orbitals
!
!
   subroutine visualize_density(this, wf)
!!
!!    Visualize density
!!    Written by Alexander C. Paul, May 2022
!!
      implicit none
!
      class(hf_visualization_task), intent(inout) :: this
!
      class(hf), intent(in) :: wf
!
      character(len=:), allocatable :: density_file_tag
!
      if (.not. this%plot_density) return
!
      density_file_tag = 'AO_density'
      call this%visualizer%plot_density(wf%ao, wf%ao_density, density_file_tag)
!
   end subroutine visualize_density
!
!
   subroutine visualize_active_density(this, wf)
!!
!!    Visualize active density
!!    Written by Ida-Marie Høyvik, Oct 2019
!!
      use visualization_class, only: visualization
!
      implicit none
!
      class(hf_visualization_task), intent(inout) :: this
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(:,:), allocatable :: density
!
      character(len=:), allocatable :: density_file_tag
!
      if (.not. this%plot_active_density) return
!
      call mem%alloc(density, wf%ao%n, wf%ao%n)
!
      call dgemm('N', 'T',                   &
                  wf%ao%n,                   &
                  wf%ao%n,                   &
                  wf%n_o,                    &
                  two,                       &
                  wf%orbital_coefficients,   &
                  wf%ao%n,                   &
                  wf%orbital_coefficients,   &
                  wf%ao%n,                   &
                  zero,                      &
                  density,                   &
                  wf%ao%n)
!
      density_file_tag = "active_" // trim(wf%name_) // '_density'
!
      call this%visualizer%plot_density(wf%ao, density, density_file_tag)
!
      call mem%dealloc(density, wf%ao%n, wf%ao%n)
!
   end subroutine visualize_active_density
!
!
end module hf_visualization_task_class
