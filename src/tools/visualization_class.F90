!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
!!    A tool to plot orbitals and densities in common file formats.
!!    This class has two routines that may be accessed outside of this module: 
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
!!       plotter = visualization(wf%ao)
!!       call plotter%initialize(wf%ao)
!!
!!    Plot something:
!!
!!       call plotter%plot_orbitals(wf%ao, orbital_coefficients, n_mo, file_tags)
!!       call plotter%plot_density(wf%ao, density, file_tags)
!!
!!       plot_orbitals  -> use it to plot orbitals by passing orbital coefficients
!!       plot_density   -> use it to plot densities (NOTE: see documentation of how 
!!                         these densities must be prepared)
!!    
!!    Deallocate arrays used by the plotter
!!
!!       call plotter%cleanup()
!!
!
   use parameters
!
   use global_out, only: output
   use global_in, only: input
   use memory_manager_class, only: mem
   use ao_tool_class, only: ao_tool
!
   type :: visualization 
!
      integer, private  :: n_ao
      integer, private  :: n_x, n_y, n_z, n_grid_points
!     The padding in x, y, and z directions around the molecule  (Angstrom)
      real(dp), private :: buffer
      real(dp), private :: dx ! Grid point spacing
      real(dp), private :: x_min, y_min, z_min
      real(dp), private :: x_max, y_max, z_max
!
      character(len=200), private :: file_format
      real(dp), dimension(:,:), allocatable, private :: aos_on_grid
!
      logical, private :: grid_in_memory
      logical, private :: overwrite_grid
!
   contains 
!
      procedure, public :: initialize                            => initialize_visualization
      procedure, public :: cleanup                               => cleanup_visualization
!
      procedure, public :: plot_orbitals => plot_orbitals_visualization
      procedure, public :: plot_density  => plot_density_visualization
!
      procedure, private :: set_up_grid                  => set_up_grid_visualization
      procedure, private :: print_grid_info              => print_grid_info_visualization
      procedure, private :: read_settings                => read_settings_visualization
      procedure, private :: write_vector                 => write_vector_visualization
      procedure, private :: write_vector_to_plt          => write_vector_to_plt_visualization
      procedure, private :: write_vector_to_cube         => write_vector_to_cube_visualization
      procedure, private :: evaluate_mos_on_grid         => evaluate_mos_on_grid_visualization
      procedure, private :: evaluate_density_on_grid     => evaluate_density_on_grid_visualization
      procedure, private :: place_grid_in_memory         => place_grid_in_memory_visualization
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
   function new_visualization(ao) result(plotter)
!!
!!    New visualization 
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Aug 2019
!!
      implicit none 
!
      type(ao_tool), intent(inout) :: ao
!
      type(visualization) :: plotter
!
!     Default settings (lengths given in Ångströms)
!
      plotter%dx              = 0.1d0 
      plotter%buffer          = 2.00d0
      plotter%file_format     = 'plt'
      plotter%grid_in_memory  = .false.
      plotter%n_ao            = ao%n
!
      call output%printf('m', ':: Visualization of orbitals and density', fs='(/t3,a)')
!
      call plotter%read_settings()
!
      call plotter%set_up_grid(ao)
!
      call plotter%print_grid_info()
!
   end function new_visualization
!
!
   subroutine initialize_visualization(plotter, ao)
!!
!!    Initialize
!!    Written by Andreas S. Skeidsvoll, Sep 2020
!!
!!    Allocates and fills an array with the AOs evaluated on the grid, if the
!!    array can be held in memory.
!!
      implicit none 
!
      class(visualization), intent(inout) :: plotter
      type(ao_tool), intent(inout) :: ao
!
      if (plotter%n_ao*plotter%n_grid_points*dp .lt. mem%get_available()/2) then
!
         plotter%grid_in_memory = .true.
         call mem%alloc(plotter%aos_on_grid, plotter%n_ao, plotter%n_grid_points)
!
         call plotter%place_grid_in_memory(ao)
!
      endif
!
   end subroutine initialize_visualization
!
!
  subroutine set_up_grid_visualization(plotter, ao)
!!
!!    Set up grid
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Jul 2019
!!
!!    Sets up the grid, that is, given the
!!    buffer and grid point spacing, it sets up a 
!!    grid around the molecule (center coordinates)
!!
      implicit none
!
      class(visualization), intent(inout) :: plotter
!
      type(ao_tool), intent(in) :: ao
!
      integer :: i
!
      real(dp), dimension(:,:), allocatable :: R_qk
!
      if (.not. plotter%overwrite_grid) then
!
         if (plotter%buffer .lt. plotter%dx) &
            call output%error_msg('in visualization tool. Buffer is smaller than grid point spacing!')
!
!        Find minimal and maximal atomic x, y and z positions in the molecule
!
         call mem%alloc(R_qk, 3, ao%get_n_centers())
!
         R_qk = ao%get_center_coordinates()
!
         plotter%x_min = R_qk(1, 1)
         plotter%y_min = R_qk(2, 1)
         plotter%z_min = R_qk(3, 1)
!
         plotter%x_max = R_qk(1, 1)
         plotter%y_max = R_qk(2, 1)
         plotter%z_max = R_qk(3, 1)
!
         do i = 2, ao%get_n_centers()
!
            plotter%x_min = min(plotter%x_min, R_qk(1, i))
            plotter%y_min = min(plotter%y_min, R_qk(2, i))
            plotter%z_min = min(plotter%z_min, R_qk(3, i))
!
            plotter%x_max = max(plotter%x_max, R_qk(1, i))
            plotter%y_max = max(plotter%y_max, R_qk(2, i))
            plotter%z_max = max(plotter%z_max, R_qk(3, i))
!
         enddo
!
         call mem%dealloc(R_qk, 3, ao%get_n_centers())
!
!        Subtract/add buffer from minimum/maximum X,Y,Z values
!
         plotter%x_min = plotter%x_min - plotter%buffer
         plotter%y_min = plotter%y_min - plotter%buffer
         plotter%z_min = plotter%z_min - plotter%buffer
!
         plotter%x_max = plotter%x_max + plotter%buffer
         plotter%y_max = plotter%y_max + plotter%buffer
         plotter%z_max = plotter%z_max + plotter%buffer
!
      endif
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
   subroutine write_vector_visualization(plotter, vector, file_name)
!!
!!    Write vector
!!    Written by Tor S. Haugland, Nov 2019
!!
!!    Wrapper for writing to different file formats.
!!
      implicit none
!
      class(visualization), intent(in) :: plotter
!
      real(dp), dimension(plotter%n_grid_points), intent(in) :: vector
      character(len=*),                           intent(in) :: file_name
!
      select case (trim(plotter%file_format))
!
         case ('plt')
            call plotter%write_vector_to_plt(vector, file_name)
!
         case ('cube')
            call plotter%write_vector_to_cube(vector, file_name)
!
         case default
            call output%error_msg("Did not recognize format '(a0)'. Use 'plt' or 'cube'.", &
                                  chars=[trim(plotter%file_format)])
!
      end select
!
   end subroutine write_vector_visualization
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
!!    Modified by Tor S. Haugaldn, Nov 2020
!!    Input vector changed from sp to dp.
!!
      use stream_file_class, only: stream_file
!
      implicit none
!
      class(visualization), intent(in) :: plotter
!
      real(dp), dimension(plotter%n_grid_points), intent(in) :: vector
!
      character(len=*), intent(in) :: file_name
!
!     Single precision
      real(sp), dimension(:), allocatable :: vector_sp
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
      n_zyx_sp(1) = int(plotter%n_z, i32)
      n_zyx_sp(2) = int(plotter%n_y, i32)
      n_zyx_sp(3) = int(plotter%n_x, i32)
!
      zyx_min_max_sp(1) = real(plotter%z_min, sp)
      zyx_min_max_sp(2) = real(plotter%z_max, sp)
      zyx_min_max_sp(3) = real(plotter%y_min, sp)
      zyx_min_max_sp(4) = real(plotter%y_max, sp)
      zyx_min_max_sp(5) = real(plotter%x_min, sp)
      zyx_min_max_sp(6) = real(plotter%x_max, sp)
!
      allocate(vector_sp(plotter%n_grid_points))
      vector_sp = real(vector, sp)
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
      call plt_file%write_(vector_sp, plotter%n_grid_points)
!
      call plt_file%close_('keep')
!
      deallocate(vector_sp)
!
   end subroutine write_vector_to_plt_visualization
!
!
   subroutine write_vector_to_cube_visualization(plotter, vector, file_name)
!!
!!    Write vector to cube file
!!    Written by Tor S. Haugland, 2019
!!
!!    Write vector of densities to a .cube file according to the specification,
!!       http://paulbourke.net/dataformats/cube/
!!
      use output_file_class, only: output_file
!
      implicit none
!
      class(visualization), intent(in) :: plotter
!
      real(dp), dimension(plotter%n_x, plotter%n_y, plotter%n_z), intent(in) :: vector
      character(len=*),                                           intent(in) :: file_name
!
      real(dp) :: x_min, y_min, z_min
      real(dp) :: x_max, y_max, z_max
      real(dp) :: dx
      real(dp) :: volume, integral
      integer:: ix, iy, iz
      type(output_file) :: cube_file
!
!     .cube files are usually in Bohr
      x_min = plotter%x_min * angstrom_to_bohr
      y_min = plotter%y_min * angstrom_to_bohr
      z_min = plotter%z_min * angstrom_to_bohr
      x_max = plotter%x_max * angstrom_to_bohr
      y_max = plotter%y_max * angstrom_to_bohr
      z_max = plotter%z_max * angstrom_to_bohr
      dx    = plotter%dx    * angstrom_to_bohr
!
!     Integrate density
      volume   = (x_max - x_min) * (y_max - y_min) * (z_max - z_min)
      integral = volume * sum(vector) / plotter%n_grid_points 
!
      cube_file = output_file(trim(file_name))
      call cube_file%open_()
!
!     Header
      call cube_file%printf('m', '! .cube file from eT (Bohr)', fs='(a)')
      call cube_file%printf('m', '! Integrated density = (f12.6)', fs='(a)', &
                           reals=[integral])
      call cube_file%printf('m', '(i12) (f12.6) (f12.6) (f12.6)', fs='(a)', &
                           ints=[0], reals=[x_min, y_min, z_min])
      call cube_file%printf('m', '(i12) (f12.6) (f12.6) (f12.6)', fs='(a)', &
                           ints=[plotter%n_x], reals=[dx, 0.0d0, 0.0d0])
      call cube_file%printf('m', '(i12) (f12.6) (f12.6) (f12.6)', fs='(a)', &
                           ints=[plotter%n_y], reals=[0.0d0, dx, 0.0d0])
      call cube_file%printf('m', '(i12) (f12.6) (f12.6) (f12.6)', fs='(a)', &
                           ints=[plotter%n_z], reals=[0.0d0, 0.0d0, dx])
!
      do ix=1, plotter%n_x
         do iy=1, plotter%n_y
            do iz=1, plotter%n_z
               call cube_file%printf('m', "(e12.6)", fs='(a,1x)', reals=[vector(ix, iy, iz)], adv=.false.)
            enddo
         enddo
         call cube_file%printf('m', "", fs='(a)')
      enddo
!
      call cube_file%close_()
!
   end subroutine write_vector_to_cube_visualization
!
!
   subroutine read_settings_visualization(plotter)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Aug 2019
!!    Modified by Tor S. Haugland Feb 2021, Add overwrite_grid
!!
      implicit none 
!
      class(visualization), intent(inout) :: plotter 
      real(dp), dimension(3) :: position
!
      call input%get_keyword('grid spacing', 'visualization', plotter%dx)
      call input%get_keyword('grid buffer', 'visualization', plotter%buffer)
      call input%get_keyword('file format', 'visualization', plotter%file_format)
!
!     Overwrite minimum/maximum X,Y,Z values in grid
!
      plotter%overwrite_grid = input%is_keyword_present('grid min', 'visualization')
      if (plotter%overwrite_grid .neqv. input%is_keyword_present('grid max', 'visualization')) &
         call output%error_msg("'grid min' and 'grid max' must be specified together.")
      
      if (plotter%overwrite_grid) then
         call input%get_array_for_keyword('grid min', 'visualization', 3, position)
         plotter%x_min = position(1)
         plotter%y_min = position(2)
         plotter%z_min = position(3)
!
         call input%get_array_for_keyword('grid max', 'visualization', 3, position)
         plotter%x_max = position(1)
         plotter%y_max = position(2)
         plotter%z_max = position(3)
      endif
!
   end subroutine read_settings_visualization
!
!
   subroutine evaluate_mos_on_grid_visualization(plotter, ao, mos_on_grid, &
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
      real(dp), dimension(plotter%n_x, plotter%n_y, &
                              plotter%n_z, n_mo), intent(out) :: mos_on_grid
!
      type(ao_tool), intent(in) :: ao
!
      real(dp), dimension(plotter%n_ao, n_mo), intent(in) :: orbital_coefficients
!
      integer :: i, j, k, mo, n_threads, thread
!
      real(dp), dimension(:,:), allocatable     :: aos_at_point
!
      real(dp)                                  :: x, y, z
!
      real(dp), external :: ddot
!
      if (plotter%grid_in_memory) then
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
                     mos_on_grid,            &
                     plotter%n_grid_points)
!
      else
!
!        In case of no omp:
!
         n_threads   = 1
         thread      = 1
!
!$       n_threads = omp_get_max_threads()
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
                  call ao%evaluate_aos_at_point(x, y, z, aos_at_point(:,thread))
!
                  do mo = 1, n_mo
!
                     mos_on_grid(i, j, k, mo) = ddot(plotter%n_ao, aos_at_point(:,thread), 1, &
                                                               orbital_coefficients(1,mo), 1)
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
   subroutine plot_orbitals_visualization(plotter, ao, orbital_coefficients, &
                                                n_mo, file_tags)
!!
!!    Plot orbitals
!!    Written by Sarai D. Folkestad and Andreas Skeidsvoll, Sep. 2019
!!
!!    Plots orbitals on grid.
!!
!!    Receives the atomic orbitals ao,
!!    the orbital coefficients for the molecular orbitals,
!!    the number of MOs to plot n_mo,
!!    and the file_tags for file names
!!
      use array_utilities, only: get_abs_max
!!
      implicit none
!
      class(visualization), intent(inout) :: plotter
!
      type(ao_tool), intent(in) :: ao
!
      integer, intent(in) :: n_mo
!
      real(dp), dimension(plotter%n_ao, n_mo), intent(in) :: orbital_coefficients
!
      character(len=200), dimension(n_mo), intent(in) :: file_tags
!
      real(dp), dimension(:,:,:,:), allocatable :: mos_on_grid
!
      character(len=200) :: file_name
!
      integer :: mo
!
      call output%printf('m', '- Plotting orbitals', fs='(/t3,a)')      
!
!     Create array of molecular orbitals at grid point
!
      call mem%alloc(mos_on_grid, plotter%n_x, plotter%n_y, plotter%n_z, n_mo)
      call plotter%evaluate_mos_on_grid(ao, mos_on_grid, orbital_coefficients, n_mo)
!
!     Write the molecular orbitals in the array to .plt files
!
      do mo = 1, n_mo
!
!        For each orbital write to file
!
         file_name = 'eT' // '.' // trim(file_tags(mo)) // '.' // trim(plotter%file_format)
!
         call plotter%write_vector(mos_on_grid(:, :, :, mo), trim(file_name))
!
      enddo
!
      call mem%dealloc(mos_on_grid, plotter%n_x, plotter%n_y, plotter%n_z, n_mo)
!
   end subroutine plot_orbitals_visualization
!
!
   subroutine plot_density_visualization(plotter, ao, density, file_tag)
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
      class(ao_tool), intent(in)           :: ao
!
      character(len=*), intent(in) :: file_tag
!
      real(dp), dimension(plotter%n_ao, plotter%n_ao), intent(in) :: density
!
      real(dp), dimension(:), allocatable :: density_on_grid_vec
!
      character(len=200) :: file_name

      call output%printf('m', '- Plotting density', fs='(/t3,a)')      
!
!     Create electron density vector
!
      call mem%alloc(density_on_grid_vec, plotter%n_grid_points)
!
      call plotter%evaluate_density_on_grid(ao, density, density_on_grid_vec)
!
!     Write vector to file
!
      file_name = 'eT' // '.' // trim(file_tag) // '.' // trim(plotter%file_format)
      call plotter%write_vector(density_on_grid_vec, trim(file_name))
!
      call mem%dealloc(density_on_grid_vec, plotter%n_grid_points)
!
   end subroutine plot_density_visualization
!
!
   subroutine evaluate_density_on_grid_visualization(plotter, ao, density, &
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
      real(dp), dimension(plotter%n_grid_points), intent(out)     :: density_on_grid_vec
!
      type(ao_tool), intent(in) :: ao
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
            density_on_grid_vec(gp) = ddot(plotter%n_ao, &
                                           plotter%aos_on_grid(1,gp), 1, &
                                           I1(1,gp), 1)
!
         enddo
!
         call mem%dealloc(I1, plotter%n_ao, plotter%n_grid_points)  
!
      else
!
!        In case of no omp:
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
                  call ao%evaluate_aos_at_point(x, y, z, aos_at_point(:,thread))
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
                  density_on_grid_vec(vector_index) = ddot(plotter%n_ao, &
                                                           aos_at_point(1,thread), 1, &
                                                           I2, 1)
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
   subroutine place_grid_in_memory_visualization(plotter, ao)
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
      type(ao_tool), intent(in)  :: ao
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
         call ao%evaluate_aos_at_point(coordinates(gp, 1), &
                                       coordinates(gp, 2), &
                                       coordinates(gp, 3), &
                                       plotter%aos_on_grid(:,gp))
!
      enddo
!$omp end parallel do
!
      call mem%dealloc(coordinates, plotter%n_grid_points, 3)
!
   end subroutine place_grid_in_memory_visualization
!
!
   subroutine cleanup_visualization(plotter)
!!
!!    Cleanup  
!!    Written by Sarai D. Folkestad, Nov 2019
!!
!!    Modified by Andreas S. Skeidsvoll, Sep 2020
!!    Changed from destructor to cleanup routine, in order to comply with
!!    Fortran standards.
!!
      implicit none 
!
      class(visualization), intent(inout) :: plotter
!
      if (plotter%grid_in_memory) call mem%dealloc(plotter%aos_on_grid, &
                                       plotter%n_ao, plotter%n_grid_points)
!
   end subroutine cleanup_visualization
!
!
end module visualization_class
