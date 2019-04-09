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
module wavefunction_class
!
!!
!!    Wavefunction class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use file_class
   use disk_manager_class
   use molecular_system_class
!
   implicit none
!
   type, abstract :: wavefunction
!
      character(len=40) :: name_
!
      real(dp) :: energy
!
      integer :: n_ao
      integer :: n_mo
      integer :: n_o
      integer :: n_v
!
      type(molecular_system), pointer :: system
!
      real(dp), dimension(:,:), allocatable :: orbital_coefficients
      real(dp), dimension(:), allocatable :: orbital_energies
!
      type(file) :: orbital_coefficients_file
      type(file) :: orbital_energies_file
!
   contains
!
      procedure :: initialize_orbital_coefficients => initialize_orbital_coefficients_wavefunction
      procedure :: initialize_orbital_energies     => initialize_orbital_energies_wavefunction
!
      procedure :: print_wavefunction_summary      => print_wavefunction_summary_wavefunction
!
      procedure :: destruct_orbital_coefficients   => destruct_orbital_coefficients_wavefunction
      procedure :: destruct_orbital_energies       => destruct_orbital_energies_wavefunction
!
      procedure :: get_ao_h_wx                     => get_ao_h_wx_wavefunction
      procedure :: get_ao_s_wx                     => get_ao_s_wx_wavefunction
      procedure :: get_ao_mu_wx                    => get_ao_mu_wx_wavefunction
      procedure :: get_ao_q_wx                     => get_ao_q_wx_wavefunction
!
      procedure :: mo_transform                    => mo_transform_wavefunction
      procedure :: mo_transform_and_save_h         => mo_transform_and_save_h_wavefunction
!
      procedure :: initialize_wavefunction_files   => initialize_wavefunction_files_wavefunction
!
      procedure :: read_orbital_coefficients                => read_orbital_coefficients_wavefunction
      procedure :: save_orbital_coefficients                => save_orbital_coefficients_wavefunction
      procedure :: read_orbital_energies                    => read_orbital_energies_wavefunction
      procedure :: save_orbital_energies                    => save_orbital_energies_wavefunction
!
      procedure :: is_restart_safe                          => is_restart_safe_wavefunction 
!
   end type wavefunction
!
contains
!
!
   subroutine initialize_orbital_coefficients_wavefunction(wf)
!!
!!    Initialize MO coefficients
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (.not. allocated(wf%orbital_coefficients)) call mem%alloc(wf%orbital_coefficients, wf%n_ao, wf%n_mo)
!
   end subroutine initialize_orbital_coefficients_wavefunction
!
!
   subroutine initialize_wavefunction_files_wavefunction(wf)
!!
!!    Initialize wavefunction files 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none 
!
      class(wavefunction) :: wf 
!
      call wf%orbital_coefficients_file%init('orbital_coefficients', 'sequential', 'unformatted')
      call wf%orbital_energies_file%init('orbital_energies', 'sequential', 'unformatted')
!
   end subroutine initialize_wavefunction_files_wavefunction
!
!
   subroutine print_wavefunction_summary_wavefunction(wf)
!!
!!    Print wavefunction summary 
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
!!    Prints information related to the wavefunction,
!!    most of which is meaningful only for a properly 
!!    converged wavefunction. Should be overwritten in 
!!    descendants if more or less or other information 
!!    is present. 
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      write(output%unit, '(/t3,a,a,a)') '- Summary of ', trim(wf%name_), ' wavefunction:'
!
!
   end subroutine print_wavefunction_summary_wavefunction
!
!
   subroutine destruct_orbital_coefficients_wavefunction(wf)
!!
!!    Destruct MO coefficients
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (allocated(wf%orbital_coefficients)) call mem%dealloc(wf%orbital_coefficients, wf%n_ao, wf%n_mo)
!
   end subroutine destruct_orbital_coefficients_wavefunction
!
!
   subroutine initialize_orbital_energies_wavefunction(wf)
!!
!!    Initialize orbital energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (.not. allocated(wf%orbital_energies)) call mem%alloc(wf%orbital_energies, wf%n_mo)
!
   end subroutine initialize_orbital_energies_wavefunction
!
!
   subroutine destruct_orbital_energies_wavefunction(wf)
!!
!!    Destruct AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (allocated(wf%orbital_energies)) call mem%dealloc(wf%orbital_energies, wf%n_mo)
!
   end subroutine destruct_orbital_energies_wavefunction
!
!
   subroutine mo_transform_and_save_h_wavefunction(wf)
!!
!!    MO transform and save h 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none 
!
      class(wavefunction) :: wf 
!
      real(dp), dimension(:,:), allocatable :: h_wx, h_pq 
!
      type(file) :: h_pq_file
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
      call mem%alloc(h_pq, wf%n_mo, wf%n_mo)
!
      call wf%get_ao_h_wx(h_wx)
      call wf%mo_transform(h_wx, h_pq)
!
      call h_pq_file%init('h_pq', 'sequential', 'unformatted')
      call disk%open_file(h_pq_file, 'write', 'rewind')
!
      write(h_pq_file%unit) h_pq 
!
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
      call mem%dealloc(h_pq, wf%n_mo, wf%n_mo)     
!
      call disk%close_file(h_pq_file)
!
   end subroutine mo_transform_and_save_h_wavefunction
!
!
   subroutine mo_transform_wavefunction(wf, X_wx, Y_pq)
!!
!!    MO transform 
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
!!    Performs MO transformation of X and saves the result in Y:
!!
!!       Y_pq = sum_wx C_wp X_wx C_xq
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)    :: X_wx 
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: Y_pq  
!
      real(dp), dimension(:,:), allocatable :: Z_wq ! = sum_x X_wx C_xq
!
      call mem%alloc(Z_wq, wf%n_ao, wf%n_mo)
!
      call dgemm('N', 'N',                 &
                  wf%n_ao,                 &
                  wf%n_mo,                 &
                  wf%n_ao,                 &
                  one,                     &
                  X_wx,                    &
                  wf%n_ao,                 &
                  wf%orbital_coefficients, & ! C_xq
                  wf%n_ao,                 &
                  zero,                    &
                  Z_wq,                    &
                  wf%n_ao)
!
      call dgemm('T', 'N',                 &
                  wf%n_mo,                 &
                  wf%n_mo,                 &
                  wf%n_ao,                 &
                  one,                     &
                  wf%orbital_coefficients, & ! C_wp 
                  wf%n_ao,                 &
                  Z_wq,                    &
                  wf%n_ao,                 &
                  zero,                    &
                  Y_pq,                    &
                  wf%n_mo)
!
      call mem%dealloc(Z_wq, wf%n_ao, wf%n_mo)
!
   end subroutine mo_transform_wavefunction
!
!
   subroutine get_ao_h_wx_wavefunction(wf, h)
!!
!!    Get AO h 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Uses the integral tool to construct the full one-electron h matrix.
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: h 
!
      type(interval) :: A_interval, B_interval
!
      integer :: x, y, A, B
!
      real(dp), dimension(:,:), allocatable :: h_AB 
!
!$omp parallel do &
!$omp private(A, B, h_AB, A_interval, B_interval, x, y) schedule(static)
      do A = 1, wf%system%n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!
            call mem%alloc(h_AB, A_interval%size, B_interval%size)
            call wf%system%construct_ao_h_wx(h_AB, A, B)
!
             do x = 1, A_interval%size
                do y = 1, B_interval%size
!
                   h(A_interval%first - 1 + x, B_interval%first - 1 + y) = h_AB(x, y)
                   h(B_interval%first - 1 + y, A_interval%first - 1 + x) = h_AB(x, y)
!
                enddo
             enddo
!
            call mem%dealloc(h_AB, A_interval%size, B_interval%size)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine get_ao_h_wx_wavefunction
!
!
   subroutine get_ao_s_wx_wavefunction(wf, s)
!!
!!    Get AO s 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Uses the integral tool to construct the full one-electron h matrix.
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: s 
!
      type(interval) :: A_interval, B_interval
!
      integer :: x, y, A, B
!
      real(dp), dimension(:,:), allocatable :: s_AB 
!
      do A = 1, wf%system%n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!
            call mem%alloc(s_AB, A_interval%size, B_interval%size)
            call wf%system%construct_ao_s_wx(s_AB, A, B)
!
             do x = 1, A_interval%size
                do y = 1, B_interval%size
!
                   s(A_interval%first - 1 + x, B_interval%first - 1 + y) = s_AB(x, y)
                   s(B_interval%first - 1 + y, A_interval%first - 1 + x) = s_AB(x, y)
!
                enddo
             enddo
!
            call mem%dealloc(s_AB, A_interval%size, B_interval%size)
!
         enddo
      enddo
!
   end subroutine get_ao_s_wx_wavefunction
!
!
   subroutine get_ao_mu_wx_wavefunction(wf, mu_X, mu_Y, mu_Z)
!!
!!    Get AO mu
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
!!    Uses the integral tool to construct the full dipole integrals
!!    for the X, Y, and Z components.
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: mu_X
      real(dp), dimension(wf%n_ao, wf%n_ao) :: mu_Y
      real(dp), dimension(wf%n_ao, wf%n_ao) :: mu_Z
!
      type(interval) :: A_interval, B_interval
!
      integer :: x, y, A, B
!
      real(dp), dimension(:,:), allocatable :: mu_AB_X 
      real(dp), dimension(:,:), allocatable :: mu_AB_Y 
      real(dp), dimension(:,:), allocatable :: mu_AB_Z 
!
      do A = 1, wf%system%n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!
            call mem%alloc(mu_AB_X, A_interval%size, B_interval%size)
            call mem%alloc(mu_AB_Y, A_interval%size, B_interval%size)
            call mem%alloc(mu_AB_Z, A_interval%size, B_interval%size)
!
            call wf%system%construct_ao_mu_wx(mu_AB_X, mu_AB_Y, mu_AB_Z, A, B)
!
             do x = 1, A_interval%size
                do y = 1, B_interval%size
!
                   mu_X(A_interval%first - 1 + x, B_interval%first - 1 + y) = mu_AB_X(x, y)
                   mu_X(B_interval%first - 1 + y, A_interval%first - 1 + x) = mu_AB_X(x, y)
!
                   mu_Y(A_interval%first - 1 + x, B_interval%first - 1 + y) = mu_AB_Y(x, y)
                   mu_Y(B_interval%first - 1 + y, A_interval%first - 1 + x) = mu_AB_Y(x, y)
!
                   mu_Z(A_interval%first - 1 + x, B_interval%first - 1 + y) = mu_AB_Z(x, y)
                   mu_Z(B_interval%first - 1 + y, A_interval%first - 1 + x) = mu_AB_Z(x, y)
!
                enddo
             enddo
!
            call mem%dealloc(mu_AB_X, A_interval%size, B_interval%size)
            call mem%dealloc(mu_AB_Y, A_interval%size, B_interval%size)
            call mem%dealloc(mu_AB_Z, A_interval%size, B_interval%size)
!
         enddo
      enddo
!
   end subroutine get_ao_mu_wx_wavefunction
!
!
   subroutine get_ao_q_wx_wavefunction(wf, q_xx, q_xy, q_xz, q_yy, q_yz, q_zz)
!!
!!    Get AO q
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
!!    Uses the integral tool to construct the quadrupole integrals.
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: q_xx
      real(dp), dimension(wf%n_ao, wf%n_ao) :: q_xy
      real(dp), dimension(wf%n_ao, wf%n_ao) :: q_xz
      real(dp), dimension(wf%n_ao, wf%n_ao) :: q_yy
      real(dp), dimension(wf%n_ao, wf%n_ao) :: q_yz
      real(dp), dimension(wf%n_ao, wf%n_ao) :: q_zz
!
      type(interval) :: A_interval, B_interval
!
      integer :: x, y, A, B
!
      real(dp), dimension(:,:), allocatable :: q_AB_xx 
      real(dp), dimension(:,:), allocatable :: q_AB_xy 
      real(dp), dimension(:,:), allocatable :: q_AB_xz 
      real(dp), dimension(:,:), allocatable :: q_AB_yy 
      real(dp), dimension(:,:), allocatable :: q_AB_yz 
      real(dp), dimension(:,:), allocatable :: q_AB_zz 
!
      do A = 1, wf%system%n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!
            call mem%alloc(q_AB_xx, A_interval%size, B_interval%size)
            call mem%alloc(q_AB_xy, A_interval%size, B_interval%size)
            call mem%alloc(q_AB_xz, A_interval%size, B_interval%size)
            call mem%alloc(q_AB_yy, A_interval%size, B_interval%size)
            call mem%alloc(q_AB_yz, A_interval%size, B_interval%size)
            call mem%alloc(q_AB_zz, A_interval%size, B_interval%size)
!
            call wf%system%ao_integrals%construct_ao_q_wx(q_AB_xx, q_AB_xy, q_AB_xz, q_AB_yy, q_AB_yz, q_AB_zz, A, B)
!
             do x = 1, A_interval%size
                do y = 1, B_interval%size
!
                   q_xx(A_interval%first - 1 + x, B_interval%first - 1 + y) = q_AB_xx(x, y)
                   q_xx(B_interval%first - 1 + y, A_interval%first - 1 + x) = q_AB_xx(x, y)
!
                   q_xy(A_interval%first - 1 + x, B_interval%first - 1 + y) = q_AB_xy(x, y)
                   q_xy(B_interval%first - 1 + y, A_interval%first - 1 + x) = q_AB_xy(x, y)
!
                   q_xz(A_interval%first - 1 + x, B_interval%first - 1 + y) = q_AB_xz(x, y)
                   q_xz(B_interval%first - 1 + y, A_interval%first - 1 + x) = q_AB_xz(x, y)
!
                   q_yy(A_interval%first - 1 + x, B_interval%first - 1 + y) = q_AB_yy(x, y)
                   q_yy(B_interval%first - 1 + y, A_interval%first - 1 + x) = q_AB_yy(x, y)
!
                   q_yz(A_interval%first - 1 + x, B_interval%first - 1 + y) = q_AB_yz(x, y)
                   q_yz(B_interval%first - 1 + y, A_interval%first - 1 + x) = q_AB_yz(x, y)
!
                   q_zz(A_interval%first - 1 + x, B_interval%first - 1 + y) = q_AB_zz(x, y)
                   q_zz(B_interval%first - 1 + y, A_interval%first - 1 + x) = q_AB_zz(x, y)
!
                enddo
             enddo
!
            call mem%dealloc(q_AB_xx, A_interval%size, B_interval%size)
            call mem%dealloc(q_AB_xy, A_interval%size, B_interval%size)
            call mem%dealloc(q_AB_xz, A_interval%size, B_interval%size)
            call mem%dealloc(q_AB_yy, A_interval%size, B_interval%size)
            call mem%dealloc(q_AB_yz, A_interval%size, B_interval%size)
            call mem%dealloc(q_AB_zz, A_interval%size, B_interval%size)
!
         enddo
      enddo
!
   end subroutine get_ao_q_wx_wavefunction
!
!
   subroutine save_orbital_coefficients_wavefunction(wf)
!!
!!    Save orbital coefficients 
!!    Written by Eirik F. Kjønstad, Oct 2018 
!!
      implicit none 
!
      class(wavefunction), intent(inout) :: wf 
!
      call disk%open_file(wf%orbital_coefficients_file, 'write', 'rewind')
!
      write(wf%orbital_coefficients_file%unit) wf%orbital_coefficients
!
      call disk%close_file(wf%orbital_coefficients_file)
!
   end subroutine save_orbital_coefficients_wavefunction
!
!
   subroutine read_orbital_coefficients_wavefunction(wf)
!!
!!    Save orbital coefficients 
!!    Written by Eirik F. Kjønstad, Oct 2018 
!!
      implicit none 
!
      class(wavefunction), intent(inout) :: wf 
!
      call wf%is_restart_safe('ground state')
!
      call disk%open_file(wf%orbital_coefficients_file, 'read', 'rewind')
!
      read(wf%orbital_coefficients_file%unit) wf%orbital_coefficients
!
      call disk%close_file(wf%orbital_coefficients_file)
!
   end subroutine read_orbital_coefficients_wavefunction
!
!
   subroutine save_orbital_energies_wavefunction(wf)
!!
!!    Save orbital energies 
!!    Written by Eirik F. Kjønstad, Oct 2018 
!!
      implicit none 
!
      class(wavefunction), intent(inout) :: wf 
!
      call disk%open_file(wf%orbital_energies_file, 'write', 'rewind')
!
      write(wf%orbital_energies_file%unit) wf%orbital_energies
!
      call disk%close_file(wf%orbital_energies_file)
!
   end subroutine save_orbital_energies_wavefunction
!
!
   subroutine read_orbital_energies_wavefunction(wf)
!!
!!    Save orbital energies 
!!    Written by Eirik F. Kjønstad, Oct 2018 
!!
      implicit none 
!
      class(wavefunction), intent(inout) :: wf 
!
      call wf%is_restart_safe('ground state')
!
      call disk%open_file(wf%orbital_energies_file, 'read', 'rewind')
!
      read(wf%orbital_energies_file%unit) wf%orbital_energies
!
      call disk%close_file(wf%orbital_energies_file)
!
   end subroutine read_orbital_energies_wavefunction
!
!
   subroutine is_restart_safe_wavefunction(wf, task)
!!
!!    Is restart safe?
!!    Written by Eirik F. Kjønstad, Mar 2019 
!!
      implicit none 
!
      class(wavefunction) :: wf 
!
      character(len=*), intent(in) :: task 
!
      call output%error_msg('Cannot restart for task ' // trim(task) // ' from abstract wavefunction ' // trim(wf%name_))
!
   end subroutine is_restart_safe_wavefunction
!
!
end module wavefunction_class
