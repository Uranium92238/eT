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
module td_cc_engine_class
!
!!
!! Time-dependent CC engine class
!! Written by Alexander C. Paul, Sarai D. Folkestad, and Eirik F. Kjønstad, 2021
!!
!
   use kinds
   use ccs_class,                         only: ccs
   use cc_engine_class,                   only: cc_engine
   use cc_task_class,                     only: cc_task
   use cc_amplitudes_task_class,          only: cc_amplitudes_task
   use cc_multipliers_task_class,         only: cc_multipliers_task
   use cc_propagation_task_class,         only: cc_propagation_task
   use global_in,                         only: input
   use memory_manager_class,              only: mem
   use eri_approximator_task_class,       only: eri_approximator_task
   use fft_task_class,                    only: fft_task
   use cc_wavefunctions_class,            only: cc_wavefunctions
!
   implicit none
!
   type, extends(cc_engine) :: td_cc_engine
!
      type(eri_approximator_task),            private :: eri_approximator
      type(cc_amplitudes_task), allocatable,  private :: ground_state_amplitudes
      type(cc_multipliers_task), allocatable, private :: ground_state_multipliers
      type(cc_propagation_task), allocatable, private :: amplitude_propagation
      type(fft_task), allocatable, private :: dipole_fft
      type(fft_task), allocatable, private :: electric_field_fft
!
!     Fast Fourier transform
      logical, private :: fft_dipole_moment, fft_electric_field
!
   contains
!
      procedure, public :: ignite => ignite_td_cc_engine
      procedure, public :: set_allowed_wfs => set_allowed_wfs_td_cc_engine
!
      procedure, private, nopass :: plot_density
!
   end type td_cc_engine
!
!
contains
!
!
   subroutine ignite_td_cc_engine(this, wf)
!!
!!    Ignite
!!    Written by Eirik F. Kjønstad, 2021
!!
      use citation_printer_class, only: eT_citations
!
      implicit none
!
      class(td_cc_engine), intent(inout) :: this
!
      class(ccs), intent(inout) :: wf
!
      call this%set_allowed_wfs()
      call this%check_wavefunctions(wf)
!
      call eT_citations%add("Time-dependent CC")
!
      this%fft_dipole_moment = input%is_keyword_present('fft dipole moment', 'cc td')
      this%fft_electric_field = input%is_keyword_present('fft electric field', 'cc td')
!
      call this%eri_approximator%execute(wf)
!
      this%ground_state_amplitudes = cc_amplitudes_task()
      call this%ground_state_amplitudes%execute(wf)
!
      this%ground_state_multipliers = cc_multipliers_task()
      call this%ground_state_multipliers%execute(wf)
!
      call wf%make_complex()
!
      this%amplitude_propagation = cc_propagation_task()
      call this%amplitude_propagation%execute(wf)
!
      if (this%fft_dipole_moment) then
!
         this%dipole_fft = fft_task('dipole moment', 'cc_propagation_dipole_moment')
         call this%dipole_fft%execute(wf)
!
      end if
!
      if (this%fft_electric_field) then
!
         this%electric_field_fft = fft_task('electric field', &
                                            'cc_propagation_electric_field')
         call this%electric_field_fft%execute(wf)
!
      end if
!
      if (input%is_keyword_present('plot cc density', 'visualization')) &
         call this%plot_density(wf)
!
      call wf%cleanup_complex()
!
   end subroutine ignite_td_cc_engine
!
!
   subroutine plot_density(wf)
!!
!!    Plot density
!!    Written by Andreas Skeidsvoll, Dec 2019
!!
!!    Reads the electron density matrices listed in cc_propagation_density_matrix_real and
!!    cc_propagation_density_matrix_imaginary, and writes the corresponding electron densities to
!!    file.
!!
!!    Based on do_visualization_gs_this by Tor S. Haugland, Nov 2019
!!
      use visualization_class, only : visualization
      use array_utilities, only: symmetric_sandwich_right_transposition
      use formatted_read_file_class, only: formatted_read_file
      use timings_class, only: timings
!
      implicit none
!
      class(ccs) :: wf
!
      type(visualization), allocatable :: plotter
!
      real(dp), dimension(:,:), allocatable :: mo_density, density
!
      type(formatted_read_file) :: density_matrix_real_file, density_matrix_imaginary_file
!
      integer :: file_count, iostat
!
      character(len=200) :: file_count_string
!
      type(timings) :: timer
!
      timer = timings('Plotting TDCC densities')
      call timer%turn_on
!
      plotter = visualization(wf%ao)
      call plotter%initialize(wf%ao)
!
      call mem%alloc(mo_density, wf%n_mo, wf%n_mo)
      call mem%alloc(density, wf%ao%n, wf%ao%n)
!
!     Plot real electron densities using density matrices on file
      density_matrix_real_file = formatted_read_file('eT.cc_propagation_density_matrix_real.out')
      call density_matrix_real_file%open_('rewind')
!
      file_count = 0
!
      do
!
         call density_matrix_real_file%read_(mo_density, wf%n_mo*wf%n_mo, io_stat=iostat)
!
         if (iostat .ne. 0) exit
!
!        D_alpha,beta = sum_pq  D_pq C_alpha,p C_beta,q
         call symmetric_sandwich_right_transposition(density,                 &
                                                     mo_density,              &
                                                     wf%orbital_coefficients, &
                                                     wf%ao%n,                 &
                                                     wf%n_mo)
!
         file_count = file_count + 1
         write(file_count_string, *) file_count
!
         call plotter%plot_density(wf%ao, density, 'cc_propagation_density_matrix_real_' &
                                                       // trim(adjustl(file_count_string)))
!
      enddo
!
      call density_matrix_real_file%close_
!
!     Plot imaginary electron densities using density matrices on file
      density_matrix_imaginary_file = formatted_read_file('eT.cc_propagation_density_matrix_imag.out')
      call density_matrix_imaginary_file%open_('rewind')
!
      file_count = 0
!
      do
!
         call density_matrix_imaginary_file%read_(mo_density, wf%n_mo*wf%n_mo, io_stat=iostat)
!
         if (iostat .ne. 0) exit
!
!        D_alpha,beta = sum_pq  D_pq C_alpha,p C_beta,q
         call symmetric_sandwich_right_transposition(density,                 &
                                                     mo_density,              &
                                                     wf%orbital_coefficients, &
                                                     wf%ao%n,                 &
                                                     wf%n_mo)
!
!        Plot density
         file_count = file_count + 1
         write(file_count_string, *) file_count
!
         call plotter%plot_density(wf%ao, density, 'cc_propagation_density_matrix_imaginary_' &
                                                       // trim(adjustl(file_count_string)))
!
      enddo
!
      call density_matrix_imaginary_file%close_
!
      call mem%dealloc(mo_density, wf%n_mo, wf%n_mo)
      call mem%dealloc(density, wf%ao%n, wf%ao%n)
!
      call plotter%cleanup()
      call timer%turn_off
!
   end subroutine plot_density
!
!
   subroutine set_allowed_wfs_td_cc_engine(this)
!!
!!    Set allowed wavefunctions
!!    Written by Sarai D. Folkestad, Feb 2022
!!
      implicit none
!
      class(td_cc_engine), intent(inout) :: this
!
      this%allowed_cc_wfs = cc_wavefunctions()
!

      call this%allowed_cc_wfs%set('ccs',  allowed=.true.)
      call this%allowed_cc_wfs%set('ccsd', allowed=.true.)
!
   end subroutine set_allowed_wfs_td_cc_engine
!
!
end module td_cc_engine_class
