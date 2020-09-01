!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
module visualization_class
!
!!
!!    Visualization class
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Aug 2019
!!
!!    A tool to plot orbitals and densities. That is, produce .plt files
!!    which may be opened in Chimera. This class has two
!!    routines that may be accessed outside of this module: 
!!
!!    Usage:
!!
!!    Include a use statement
!!
!!       use visualization_class, only : visualization
!!
!!    Declare and initialize a tool
!!
!!       type(visualization) :: plotter
!!
!!       ...
!!
!!       plotter = visualization(wf%system)
!!
!!    Plot something:
!!
!!       call plotter%plot_orbitals(wf%system, orbital_coefficients, wf%n_ao, n_mo, file_tags)
!!       call plotter%plot_density(wf%system, wf%n_ao, D, file_tag)
!!
!!       plot_orbitals  -> use it to plot orbitals by passing orbital coefficients
!!       plot_density   -> use it to plot densities (NOTE: see documentation of how 
!!                         these densities must be prepared)
!!    
!!
!
   use parameters
!
   use molecular_system_class, only: molecular_system
   use global_out, only: output
   use global_in, only: input
   use memory_manager_class, only: mem
!
   type :: visualization 
!
      integer  :: n_ao
!
      real(dp) :: buffer ! The padding in x, y, and z directions around the molecule  (Angstrom)
!
      integer  :: n_x, n_y, n_z, n_grid_points
      real(dp) :: dx ! Grid point spacing
      real(dp) :: x_min, y_min, z_min
      real(dp) :: x_max, y_max, z_max
!
      real(dp), dimension(:,:), allocatable :: aos_on_grid
!
      logical :: grid_in_memory
!
   contains 
!
      procedure :: plot_orbitals => plot_orbitals_visualization
      procedure :: plot_density  => plot_density_visualization
!
      procedure, private :: set_up_grid                  => set_up_grid_visualization
      procedure, private :: print_grid_info              => print_grid_info_visualization
      procedure, private :: read_settings                => read_settings_visualization
      procedure, private :: write_vector_to_plt          => write_vector_to_plt_visualization
      procedure, private :: evaluate_mos_on_grid         => evaluate_mos_on_grid_visualization
      procedure, private :: evaluate_density_on_grid     => evaluate_density_on_grid_visualization
      procedure, private :: place_grid_in_memory         => place_grid_in_memory_visualization
!
      final :: destructor_visualization
!
   end type visualization
!
!
   interface visualization
!
      procedure :: new_visualization
!
   end interface visualization 
!
!
contains
!
!
   function new_visualization(system, n_ao) result(plotter)
!!
!!    New visualization 
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Aug 2019
!!
      implicit none 
!
      type(molecular_system), intent(inout) :: system
      type(visualization) :: plotter
!
      integer :: n_ao
!
!     Default settings (lengths given in Ångströms)
!
      plotter%dx              = 0.1d0 
      plotter%buffer          = 2.00d0
      plotter%grid_in_memory  = .false.
      plotter%n_ao            = n_ao
!
      call output%printf('m', ':: Visualization of orbitals and density', fs='(/t3,a)')
!
      call plotter%read_settings()
!
!     Set up grid
!
      call plotter%set_up_grid(system)
!
      call plotter%print_grid_info()
!
      call system%set_basis_info()
!
!     Evaluate aos on grid if we can keep it in memory
!
     if (plotter%n_ao*plotter%n_grid_points*dp .lt. mem%get_available()/2) then
!
        plotter%grid_in_memory = .true.
        call mem%alloc(plotter%aos_on_grid, plotter%n_ao, plotter%n_grid_points)
!
        call plotter%place_grid_in_memory(system)
!
     endif
!
   end function new_visualization
!
!
  subroutine set_up_grid_visualization(plotter, system)
!!
!!    Set up grid
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Jul 2019
!!
!!    Sets up the grid, that is, given the
!!    buffer and grid point spacing, it sets up a 
!!    grid around the molecule (system)
!!
!!
      implicit none
!
      class(visualization), intent(inout) :: plotter
      type(molecular_system), intent(in)      :: system
!
      integer :: i
!
      if (plotter%buffer .lt. plotter%dx) &
         call output%error_msg('in visualization tool. Buffer is smaller than grid point spacing!')
!
!     Find minimal and maximal atomic x, y and z positions in the molecule
!
      plotter%x_min = system%atoms(1)%x
      plotter%y_min = system%atoms(1)%y
      plotter%z_min = system%atoms(1)%z
!
      plotter%x_max = system%atoms(1)%x
      plotter%y_max = system%atoms(1)%y
      plotter%z_max = system%atoms(1)%z
!
      do i = 2, system%n_atoms
!
         plotter%x_min = min(plotter%x_min, system%atoms(i)%x)
         plotter%y_min = min(plotter%y_min, system%atoms(i)%y)
         plotter%z_min = min(plotter%z_min, system%atoms(i)%z)
!
         plotter%x_max = max(plotter%x_max, system%atoms(i)%x)
         plotter%y_max = max(plotter%y_max, system%atoms(i)%y)
         plotter%z_max = max(plotter%z_max, system%atoms(i)%z)
!
      enddo
!
!     Subtract/add buffer from minimum/maximum X,Y,Z values
!
      plotter%x_min = plotter%x_min - plotter%buffer
      plotter%y_min = plotter%y_min - plotter%buffer
      plotter%z_min = plotter%z_min - plotter%buffer
!
      plotter%x_max = plotter%x_max + plotter%buffer
      plotter%y_max = plotter%y_max + plotter%buffer
      plotter%z_max = plotter%z_max + plotter%buffer
!
!     Determine the number of grid points in x, y, and z directions
!
      plotter%n_x = nint((plotter%x_max - plotter%x_min)/plotter%dx)
      plotter%n_y = nint((plotter%y_max - plotter%y_min)/plotter%dx)
      plotter%n_z = nint((plotter%z_max - plotter%z_min)/plotter%dx)
!
!     Update x_max, y_max, z_max
!
      plotter%x_max = plotter%x_min + (plotter%n_x - 1)*plotter%dx
      plotter%y_max = plotter%y_min + (plotter%n_y - 1)*plotter%dx
      plotter%z_max = plotter%z_min + (plotter%n_z - 1)*plotter%dx
! 
      plotter%n_grid_points = (plotter%n_x)*(plotter%n_y)*(plotter%n_z)
!
   end subroutine set_up_grid_visualization
!
!
   subroutine print_grid_info_visualization(plotter)
!!
!!    Print grid information
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Aug 2019
!!
      implicit none
!
      class(visualization), intent(in) :: plotter
!
      call output%printf('n', 'Grid information              x             y    &
                         &         z       ', fs='(/t6, a)')
!
      call output%print_separator(pl='normal', fs='(t6, a)', n=66, symbol='-')
!
      call output%printf('n', 'First (A):                (f8.2)      (f8.2)     &
                         & (f8.2)', reals=[plotter%x_min, plotter%y_min, &
                         plotter%z_min], fs='(t6, a)')
!
      call output%printf('n', 'Last (A):                 (f8.2)      (f8.2)     &
                         & (f8.2)', &
                         reals=[plotter%x_min + (plotter%n_x-1)*plotter%dx, &
                         plotter%y_min + (plotter%n_y-1)*plotter%dx, &
                         plotter%z_min + (plotter%n_z-1)*plotter%dx], fs='(t6, a)')
!
      call output%printf('n', 'Number of grid points:  (i8)      (i8)       (i8)', &
                         ints=[plotter%n_x, plotter%n_y, plotter%n_z], fs='(t6, a)')
!
      call output%print_separator(pl='normal', fs='(t6, a)', n=66, symbol='-')
!
   end subroutine print_grid_info_visualization
!  
!
   subroutine write_vector_to_plt_visualization(plotter, vector, file_name)
!!
!!    Write vector to plt file
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Aug 2019
!!
!!    Modified by Alexander C. Paul, May 2020
!!    Handles the chimera format and writes to file
!!
!
      use stream_file_class, only: stream_file
!
      implicit none
!
      class(visualization), intent(in) :: plotter
!
      real(sp), dimension(plotter%n_grid_points), intent(in) :: vector
!
      character(len=*), intent(in) :: file_name
!
!     Single precision
      integer(i32), dimension(2) :: chimera_ints
      integer(i32), dimension(3) :: n_zyx_sp
      real(sp), dimension(6)     :: zyx_min_max_sp
!
      type(stream_file) :: plt_file
!
      chimera_ints(1) = 3   ! Required integer that must always be 3
      chimera_ints(2) = 200 ! Required integer that can be anything
!
!     Type conversion to single precision for the Grid dimensions and points
!
      n_zyx_sp(1) = int(plotter%n_z,i32)
      n_zyx_sp(2) = int(plotter%n_y,i32)
      n_zyx_sp(3) = int(plotter%n_x,i32)
!
      zyx_min_max_sp(1) = real(plotter%z_min,4)
      zyx_min_max_sp(2) = real(plotter%z_max,4)
      zyx_min_max_sp(3) = real(plotter%y_min,4)
      zyx_min_max_sp(4) = real(plotter%y_max,4)
      zyx_min_max_sp(5) = real(plotter%x_min,4)
      zyx_min_max_sp(6) = real(plotter%x_max,4)
!
      plt_file = stream_file(trim(file_name))
!
      call plt_file%open_('write', 'rewind')
!
!     Write in the chimera format
!
      call plt_file%write_(chimera_ints, 2)
      call plt_file%write_(n_zyx_sp, 3)
      call plt_file%write_(zyx_min_max_sp, 6)
      call plt_file%write_(vector, size(vector))
!
      call plt_file%close_('keep')
!
   end subroutine write_vector_to_plt_visualization
!
!
   subroutine read_settings_visualization(plotter)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Aug 2019
!!
      implicit none 
!
      class(visualization), intent(inout) :: plotter 
!
      call input%get_keyword_in_section('grid spacing', 'visualization', plotter%dx)
      call input%get_keyword_in_section('grid buffer', 'visualization', plotter%buffer)
!
   end subroutine read_settings_visualization
!
!
   subroutine evaluate_mos_on_grid_visualization(plotter, system, mos_on_grid, &
                                                         orbital_coefficients, n_mo)
!!
!!    Evaluate MOs on grid
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Jul 2019
!!
!!    Creates the values at each grid point of the n_mo orbitals
!!    given by orbital_coefficients.
!!    This is done by contracting the vector of AO values on the grid
!!    with the orbital coefficients of the MOs.   
!!
!!    This routine may be used for canonical orbitals 
!!    NTOs, CNTOs and so on.
!!
!
      use omp_lib
!
      implicit none
!
      class(visualization), intent(in)  :: plotter
!
      integer, intent(in) :: n_mo
!
      real(sp), dimension(plotter%n_x, plotter%n_y, &
                              plotter%n_z, n_mo), intent(out) :: mos_on_grid
!
      type(molecular_system), intent(in) :: system
!
      real(dp), dimension(plotter%n_ao, n_mo), intent(in) :: orbital_coefficients
!
      integer :: i, j, k, mo, n_threads, thread
!
      real(dp), dimension(:,:), allocatable     :: aos_at_point
      real(dp), dimension(:,:,:,:), allocatable :: mos_on_grid_dp
!
      real(dp)                                  :: x, y, z
!
      real(dp), external :: ddot
!
      if (plotter%grid_in_memory) then
!
         call mem%alloc(mos_on_grid_dp, plotter%n_x, plotter%n_y, plotter%n_z, n_mo)
!
         call dgemm('T', 'N',                & 
                     plotter%n_grid_points,  &
                     n_mo,                   &
                     plotter%n_ao,           &
                     one,                    &
                     plotter%aos_on_grid,    &
                     plotter%n_ao,           &
                     orbital_coefficients,   &
                     plotter%n_ao,           &
                     zero,                   &
                     mos_on_grid_dp,         &
                     plotter%n_grid_points)
!
!$omp parallel do private (mo, k, j, i)
         do mo = 1, n_mo
            do k = 1, plotter%n_z
               do j = 1, plotter%n_y
                  do i = 1, plotter%n_x
!
                     mos_on_grid(i, j, k, mo) = real(mos_on_grid_dp(i, j, k, mo), kind=sp)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(mos_on_grid_dp, plotter%n_x, plotter%n_y, plotter%n_z, n_mo)
!
      else
!
!     In case of no omp:
!
      n_threads   = 1
      thread      = 1
!
!$    n_threads = omp_get_max_threads()
!
         call mem%alloc(aos_at_point, plotter%n_ao, n_threads)
!
!$omp parallel do private(k, j, i, x, y, z, mo, thread)
         do k = 1, plotter%n_z
            do j = 1, plotter%n_y
               do i = 1, plotter%n_x
!
!$                thread = omp_get_thread_num() + 1
!
                  x = plotter%x_min + real((i-1), dp)*plotter%dx
                  y = plotter%y_min + real((j-1), dp)*plotter%dx
                  z = plotter%z_min + real((k-1), dp)*plotter%dx
!
                  call system%evaluate_aos_at_point(x, y, z, aos_at_point(:,thread), plotter%n_ao)
!
                  do mo = 1, n_mo
!
                     mos_on_grid(i, j, k, mo) = real(ddot(plotter%n_ao, aos_at_point(:,thread), &
                         1, orbital_coefficients(1,mo), 1), kind=sp)
!
                  enddo
!
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(aos_at_point, plotter%n_ao, n_threads)
!
      endif
!
   end subroutine evaluate_mos_on_grid_visualization
!
!
   subroutine plot_orbitals_visualization(plotter, system, orbital_coefficients, &
                                                n_mo, file_tags)
!!
!!    Plot orbitals
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Sep. 2019
!!
!!    Plots orbitals on grid.
!!
!!    Receives the number of MOs to plot n_mo, 
!!    the orbital coefficients for these orbitals,
!!    the molecular system, and the file_tags for
!!    file names
!!
      use array_utilities, only: get_abs_max
!!
      implicit none
!
      class(visualization), intent(inout) :: plotter
!
      type(molecular_system), intent(in) :: system
!
      integer, intent(in) :: n_mo
!
      real(dp), dimension(plotter%n_ao, n_mo), intent(in) :: orbital_coefficients
!
      character(len=200), dimension(n_mo), intent(in) :: file_tags
!
      real(sp), dimension(:,:,:,:), allocatable :: mos_on_grid
!
      character(len=200) :: file_name
!
      integer :: mo
!
      call output%printf('m', '- Plotting orbitals', fs='(/t3,a)')      
!
!     Create array of molecular orbitals at grid point
!
      allocate(mos_on_grid(plotter%n_x, plotter%n_y, plotter%n_z, n_mo))
      call plotter%evaluate_mos_on_grid(system, mos_on_grid, orbital_coefficients, n_mo)
!
!     Write the molecular orbitals in the array to .plt files
!
      do mo = 1, n_mo
!
!        For each orbital write to file
!
         file_name = 'eT' // '.' // trim(file_tags(mo)) // '.plt'
!
         call plotter%write_vector_to_plt(mos_on_grid(:, :, :, mo), trim(file_name))
!
      enddo
!
      deallocate(mos_on_grid)
!
   end subroutine plot_orbitals_visualization
!
!
   subroutine plot_density_visualization(plotter, system, density, file_tag)
!!
!!    Plot density 
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Sep 2019
!!
!!    Plot density (e.g. AO density,  CC densities or density
!!    differences)
!!
!!    NOTE:
!!
!!    The density passed to this routine has dimension
!!     (n_ao x n_ao) and
!!
!!    For HF,  D^AO is passed to the routine.
!!
!!    For coupled cluster densities we have
!!
!!       D_alpha,beta = (sum_pq  D_pq C_alpha,p C_beta,q) 
!! 
!!    is passed to the routine.
!!
      implicit none 
!
      class(visualization), intent(in)     :: plotter
      class(molecular_system), intent(in)  :: system
!
      character(len=*), intent(in) :: file_tag
!
      real(dp), dimension(plotter%n_ao, plotter%n_ao), intent(in) :: density
!
      real(sp), dimension(:), allocatable :: density_on_grid_vec
!
      character(len=200) :: file_name

      call output%printf('m', '- Plotting density', fs='(/t3,a)')      
!
!     Create electron density vector
!
      allocate(density_on_grid_vec(plotter%n_grid_points))
!
      call plotter%evaluate_density_on_grid(system, density, density_on_grid_vec)
!
!     Write vector to .plt file
!
      file_name = 'eT' // '.' // trim(file_tag) // '.plt'
      call plotter%write_vector_to_plt(density_on_grid_vec, trim(file_name))
!
      deallocate(density_on_grid_vec)
!
   end subroutine plot_density_visualization
!
!
   subroutine evaluate_density_on_grid_visualization(plotter, system, density, &
                                                      density_on_grid_vec)
!!
!!    Evaluate density on grid
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Aug 2019
!!
!!    Calculates the expectation value of the density for each grid point:
!!
!!      rho(r) = sum_pq D_pq phi_p(r) phi_q(r) 
!!           
!!    See eqn. (2.7.33) in Molecular Electronic Structure Theory
!!
!!    For the HF density (D_pq = delta_pq nu_q), the expression reduces to 
!!
!!       sum_alpha,beta xi_alpha(r) D^AO_alpha,beta xi_beta(r)
!!
!!    and the AO density must be passed to this routine.
!!
!!    For coupled cluster densities we have
!!     
!!       sum_pq D_pq phi_p(r) phi_q(r) 
!!       = sum_alpha,beta (sum_pq  D_pq C_alpha,p C_beta,q) xi_beta(r) xi_alpha(r)
!!       = sum_alpha,beta D_alpha,beta xi_beta(r) xi_alpha(r)
!!
!!    and 
!!
!!       D_alpha,beta = (sum_pq  D_pq C_alpha,p C_beta,q) 
!! 
!!    is passed to the routine.
!
      use omp_lib
!
      implicit none
!
      class(visualization), intent(in)  :: plotter
!
      real(dp), dimension(plotter%n_ao, plotter%n_ao), intent(in) :: density
      real(sp), dimension(plotter%n_grid_points), intent(out)     :: density_on_grid_vec
!
      type(molecular_system), intent(in) :: system
!
      integer                                :: i, j, k, vector_index, gp, n_threads, thread
      real(dp), dimension(:), allocatable    :: I2
      real(dp), dimension(:,:), allocatable  :: aos_at_point, I1
      real(dp)                               :: x, y, z
!
      real(dp), external :: ddot
!
      if (plotter%grid_in_memory) then
!     
         call mem%alloc(I1, plotter%n_ao, plotter%n_grid_points)
!
         call dgemm('N','N',              &
                  plotter%n_ao,           &
                  plotter%n_grid_points,  &
                  plotter%n_ao,           & 
                  one,                    &
                  density,                &
                  plotter%n_ao,           & 
                  plotter%aos_on_grid,    &
                  plotter%n_ao,           & 
                  zero,                   &
                  I1,                     &
                  plotter%n_ao)   
!
         do gp = 1, plotter%n_grid_points 
!
            density_on_grid_vec(gp) = &
                     real(ddot(plotter%n_ao, plotter%aos_on_grid(1,gp), 1, I1(1,gp), 1), &
                          kind=sp)
!
         enddo
!
         call mem%dealloc(I1, plotter%n_ao, plotter%n_grid_points)  
!
      else
!
!     In case of no omp:
!
      n_threads   = 1
      thread      = 1
!
!$    n_threads = omp_get_max_threads()
!
         call mem%alloc(aos_at_point, plotter%n_ao, n_threads)
         call mem%alloc(I2, plotter%n_ao)
!
!$omp parallel do private(k, j, i, x, y, z, I2, vector_index, thread)
         do k = 1, plotter%n_z
            do j = 1, plotter%n_y
               do i = 1, plotter%n_x
!
!$                thread = omp_get_thread_num() + 1
!
                  x = plotter%x_min + real((i-1), dp)*plotter%dx
                  y = plotter%y_min + real((j-1), dp)*plotter%dx
                  z = plotter%z_min + real((k-1), dp)*plotter%dx
!
                  call system%evaluate_aos_at_point(x, y, z, aos_at_point(:,thread), plotter%n_ao)
!
                  call dgemv('N',                        &
                             plotter%n_ao,               &
                             plotter%n_ao,               & 
                             one,                        &
                             density,                    &
                             plotter%n_ao,               & 
                             aos_at_point(1,thread),     &
                             1,                          &
                             zero,                       &
                             I2,                         &
                             1)
!
                  vector_index = (k-1)*plotter%n_y*plotter%n_x + (j-1)*plotter%n_x + i
                  density_on_grid_vec(vector_index) = real(ddot(plotter%n_ao, &
                                                    aos_at_point(1,thread), 1, I2, 1), kind=sp)
!
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(aos_at_point, plotter%n_ao, n_threads)
         call mem%dealloc(I2, plotter%n_ao)
!
      endif
!
   end subroutine evaluate_density_on_grid_visualization
!
!
   subroutine place_grid_in_memory_visualization(plotter, system)
!!
!!    Place grid in memory
!!    Written by Sarai D. Folkestad, Nov 2019
!!
!!    Places the aos evaluated on the grid in memory.
!!
!!    The routine is only called when the array 
!!    can be placed in memory with a good marigin.
!!
      implicit none
!
      class(visualization), intent(inout) :: plotter
!
      type(molecular_system), intent(in)  :: system
!
      real(dp), dimension(:,:), allocatable :: coordinates
!
      integer :: i, j, k, gp
!
      call output%printf('m', '- Placing the AOs evaluated on the grid in memory', &
                         fs='(/t3, a)')
!
      call mem%alloc(coordinates, plotter%n_grid_points, 3) ! x, y, z
!
!$omp parallel do private(k, j, i, gp)
      do k = 1, plotter%n_z
         do j = 1, plotter%n_y
            do i = 1, plotter%n_x
!
               gp = (k-1)*plotter%n_y*plotter%n_x + (j-1)*plotter%n_x + i
!
               coordinates(gp, 1) = plotter%x_min + real((i-1), dp)*plotter%dx
               coordinates(gp, 2) = plotter%y_min + real((j-1), dp)*plotter%dx
               coordinates(gp, 3) = plotter%z_min + real((k-1), dp)*plotter%dx
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do private(gp)
      do gp = 1, plotter%n_grid_points
!
         call system%evaluate_aos_at_point(coordinates(gp, 1), &
                                             coordinates(gp, 2), &
                                             coordinates(gp, 3), &
                                             plotter%aos_on_grid(:,gp), &
                                             plotter%n_ao)
!
      enddo
!$omp end parallel do
!
      call mem%dealloc(coordinates, plotter%n_grid_points, 3)
!
   end subroutine place_grid_in_memory_visualization
!
!
   subroutine destructor_visualization(plotter)
!!
!!    Destructor  
!!    Written by Sarai D. Folkestad, Nov 2019 
!!
      implicit none 
!
      type(visualization) :: plotter
!
      if (plotter%grid_in_memory) call mem%dealloc(plotter%aos_on_grid, &
                                       plotter%n_ao, plotter%n_grid_points)
!
   end subroutine destructor_visualization
!
!
end module visualization_class
