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
module cuhf_class
!
!!
!!    Constrained unrestricted Hartree-Fock (CUHF) class module
!!    Written by Sarai D. Folkestad 2020
!!
!!    See Tsuchimochi and Scuseria, J. Chem. Phys. 133, 141102 (2010)
!!
!
   use parameters
!
   use uhf_class,            only: uhf
   use memory_manager_class, only: mem
!
   implicit none
!
   type, extends(uhf) :: cuhf
!
      real(dp), dimension(:,:), allocatable :: ao_F_cs
!
   contains
!
      procedure :: construct_closed_shell_fock        => construct_closed_shell_fock_cuhf
      procedure :: construct_spin_fock_contribution   => construct_spin_fock_contribution_cuhf
      procedure :: update_fock_and_energy             => update_fock_and_energy_cuhf
      procedure :: cleanup                            => cleanup_cuhf
      procedure :: initialize_fock                    => initialize_fock_cuhf
      procedure :: project_out_cv_and_vc              => project_out_cv_and_vc_cuhf
!
   end type cuhf
!
   interface cuhf
!
      procedure :: new_cuhf
!
   end interface cuhf
!
contains
!
!
   function new_cuhf() result(wf)
!!
!!    New CUHF
!!    Written by Sarai D. Folkestad 2020
!!
      implicit none
!
      type(cuhf) :: wf
!
      wf%name_ = 'cuhf'
      wf%cumulative_fock = .false.
      wf%fractional_uniform_valence = .false.
!
      call wf%read_settings()
!
      call wf%print_banner()
!
   end function new_cuhf
!
!
   subroutine construct_closed_shell_fock_cuhf(wf, cumulative)
!!
!!    Construct closed shell fock
!!    Written by Sarai D. Folkestad, 2020
!!
!!    Constructs
!!
!!       F^cs_wx = h_wx + (g_wxyz - 1/2 g_wzyx) D_zy
!!
      implicit none
!
      class(cuhf)          :: wf
!
      real(dp), dimension(:,:), allocatable :: G, D_diff
!
      logical :: cumulative
!
      if (cumulative) then
!
         call mem%alloc(D_diff, wf%ao%n, wf%ao%n)
         call mem%alloc(G, wf%ao%n, wf%ao%n)
!
         call dcopy(wf%ao%n**2, wf%ao_density, 1, D_diff, 1)
!
         call daxpy(wf%ao%n**2, -one, wf%previous_ao_density(1,1,1), 1, D_diff, 1)
         call daxpy(wf%ao%n**2, -one, wf%previous_ao_density(1,1,2), 1, D_diff, 1)
!
         call wf%construct_ao_G(D_diff, G)
         call daxpy(wf%ao%n**2, half, G, 1, wf%ao_F_cs, 1)
!
         call mem%dealloc(D_diff, wf%ao%n, wf%ao%n)
         call mem%dealloc(G, wf%ao%n, wf%ao%n)
!
      else
!
         call wf%construct_ao_G(wf%ao_density, wf%ao_F_cs)
         call dscal(wf%ao%n**2, half, wf%ao_F_cs, 1)
         call daxpy(wf%ao%n**2, one, wf%ao%h, 1, wf%ao_F_cs, 1)
!
      endif
!
   end subroutine construct_closed_shell_fock_cuhf
!
!
   subroutine construct_spin_fock_contribution_cuhf(wf)
!!
!!    Construct spin fock contribution
!!    Written by Sarai D. Folkestad, 2020
!!
!!    Constructs
!!
!!       F^c_wz = 1/2 g_wxyz(D^alpha_xy - D^beta_xy)
!!
!!    Constributions are stored in alpha/beta fock matrices
!!
!!       F^alpha_wz += F^c_wz
!!       F^beta_wz  -= F^c_wz
!!
!
      use omp_lib
      use array_utilities, only: copy_and_scale, get_abs_max, zero_array
      use reordering, only: symmetric_sum
!
      implicit none
!
      class(cuhf) :: wf
!
      real(dp), dimension(:,:), allocatable :: D_spin_diff
!
      integer :: thread = 0, n_threads = 1
!
      real(dp), dimension(:,:), allocatable :: F, shp_density_schwarz
!
      real(dp) :: exchange_thr, precision_thr
!
      real(dp) :: max_D_schwarz, max_eri_schwarz
!
      exchange_thr  = wf%exchange_threshold
      precision_thr = wf%ao%get_libint_epsilon()
!
      call mem%alloc(D_spin_diff, wf%ao%n, wf%ao%n)
!
      call dcopy(wf%ao%n**2, wf%ao_density_a, 1, D_spin_diff, 1)
      call daxpy(wf%ao%n**2, -one, wf%ao_density_b, 1, D_spin_diff, 1)
!
      call zero_array(wf%ao_fock_b, wf%ao%n**2)
      call zero_array(wf%ao_fock_a, wf%ao%n**2)
!
!$    n_threads = omp_get_max_threads()
!
      call mem%alloc(F, wf%ao%n, wf%ao%n*n_threads) ! [F(thread 1) F(thread 2) ...]
      call zero_array(F, (wf%ao%n**2)*n_threads)
!
      max_eri_schwarz = get_abs_max(wf%ao%cs_eri_max, wf%ao%n_sh*(wf%ao%n_sh + 1)/2)
!
      call mem%alloc(shp_density_schwarz, wf%ao%n_sh, wf%ao%n_sh)
      call wf%construct_shp_density_schwarz(shp_density_schwarz, D_spin_diff)
      max_D_schwarz = get_abs_max(shp_density_schwarz, wf%ao%n_sh**2)
!
      call wf%construct_exchange_ao_G(F, D_spin_diff, n_threads, max_D_schwarz,                  &
                                      max_eri_schwarz, shp_density_schwarz, wf%ao%n_sig_eri_shp, &
                                      exchange_thr, precision_thr, wf%ao%shells) ! -1/2 D_wx g_yxwz
!
      call mem%dealloc(shp_density_schwarz, wf%ao%n_sh, wf%ao%n_sh)
      call mem%dealloc(D_spin_diff, wf%ao%n, wf%ao%n)
!
      do thread = 1, n_threads
!
         call daxpy(wf%ao%n**2, half, F(1, (thread - 1)*wf%ao%n + 1), 1, wf%ao_fock_a, 1)
!
      enddo
!
      call mem%dealloc(F, wf%ao%n, wf%ao%n*n_threads)
!
      call symmetric_sum(wf%ao_fock_a, wf%ao%n)
      call dscal(wf%ao%n**2, half, wf%ao_fock_a, 1)
!
      call copy_and_scale(-one, wf%ao_fock_a, wf%ao_fock_b, wf%ao%n**2)
!
   end subroutine construct_spin_fock_contribution_cuhf
!
!
   subroutine update_fock_and_energy_cuhf(wf, cumulative)
!!
!!    Update Fock and energy
!!    Written by Sarai D. Folkestad, 2021
!!
!
      implicit none
!
      class(cuhf), intent(inout) :: wf
      logical, intent(in)  :: cumulative
!
      call wf%construct_closed_shell_fock(cumulative)
      call wf%construct_spin_fock_contribution()
!
      call wf%project_out_cv_and_vc(wf%ao_fock_a)
      call wf%project_out_cv_and_vc(wf%ao_fock_b)
!
      call daxpy(wf%ao%n**2, one, wf%ao_F_cs, 1, wf%ao_fock_b, 1)
      call daxpy(wf%ao%n**2, one, wf%ao_F_cs, 1, wf%ao_fock_a, 1)
!
      call wf%calculate_uhf_energy(wf%ao%h)
!
   end subroutine update_fock_and_energy_cuhf
!
!
   subroutine cleanup_cuhf(wf)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(cuhf) :: wf
!
      call wf%save_ao_density()
!
      call wf%destruct_orbital_energies()
      call wf%destruct_orbital_coefficients()
      call wf%destruct_fock()
      call wf%destruct_ao_density()
!
      call mem%dealloc(wf%ao_F_cs, wf%ao%n, wf%ao%n)
!
      call wf%destruct_mo_fock()
!
      deallocate(wf%ao)
!
   end subroutine cleanup_cuhf
!
!
   subroutine initialize_fock_cuhf(wf)
!!
!!    Initialize Fock
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the AO Fock matrix (or matrices).
!!    In spin-unrestricted wavefunctions, this alpha and beta matrices,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(cuhf) :: wf
!
      call wf%initialize_ao_fock()
!
      call wf%initialize_ao_fock_a()
      call wf%initialize_ao_fock_b()
!
      call mem%alloc(wf%ao_F_cs, wf%ao%n, wf%ao%n)
!
   end subroutine initialize_fock_cuhf
!
!
   subroutine project_out_cv_and_vc_cuhf(wf, F_correction)
!!
!!    Project out cv and vc
!!    Written by Sarai D. Folkestad, 2020
!!
!!    Enforces the constraint to remove spin contamination.
!!    Projects out the core-virtual and virtual-core blocks
!!    of F^c
!!
!!    See eqs (20a) - (20c) of J. Chem. Phys. 133, 141102 (2010), https://doi.org/10.1063/1.3503173
!!
!
      use array_utilities, only: symmetric_sandwich, symmetric_sandwich_right_transposition
!
      implicit none
!
      class(cuhf) :: wf
      real(dp), dimension(wf%ao%n, wf%ao%n) :: F_correction
!
      real(dp), dimension(:), allocatable :: occupation_numbers
      real(dp), dimension(:,:), allocatable :: C_natural_occupation_basis
      real(dp), dimension(:,:), allocatable :: F_natural_occupation_basis
      real(dp), dimension(:,:), allocatable :: SC
      integer :: i, n_open, first_open, j
!
      call mem%alloc(occupation_numbers, wf%n_mo)
      call mem%alloc(C_natural_occupation_basis, wf%ao%n, wf%n_mo)
!
      call wf%construct_natural_occupation_number_basis(C_natural_occupation_basis, occupation_numbers)
!
      first_open  = -1
      n_open      = 0
!
      do i = 1, wf%n_mo
!
         if (abs(occupation_numbers(i) - half) .lt. 1.0d-8) then
            if (first_open == -1) first_open = i
            n_open = n_open + 1
         endif
!
      enddo
!
!     Transform to the natural occupation number basis
      call mem%alloc(F_natural_occupation_basis, wf%n_mo, wf%n_mo)
      call symmetric_sandwich(F_natural_occupation_basis, &
                              F_correction,               &
                              C_natural_occupation_basis, &
                              wf%ao%n, wf%n_mo)
!
!
!     Project out components
      if (first_open > zero) then
         do i = 1, first_open - 1
            do j = first_open + n_open, wf%n_mo
!
               F_natural_occupation_basis(i,j)  = zero
               F_natural_occupation_basis(j, i) = zero
!
            enddo
         enddo
      endif
!
!
      call mem%alloc(SC, wf%ao%n, wf%n_mo)
!
      call dgemm('n', 'n', &
                  wf%ao%n, &
                  wf%n_mo, &
                  wf%ao%n, &
                  one,     &
                  wf%ao%s, &
                  wf%ao%n, &
                  C_natural_occupation_basis, &
                  wf%ao%n, &
                  zero,    &
                  SC,      &
                  wf%ao%n)
!
!     Transform back to the AO basis
      call symmetric_sandwich_right_transposition(F_correction, &
                                                  F_natural_occupation_basis, &
                                                  SC, wf%n_mo, wf%ao%n)
!
      call mem%dealloc(C_natural_occupation_basis, wf%ao%n, wf%n_mo)
      call mem%dealloc(F_natural_occupation_basis, wf%n_mo, wf%n_mo)
      call mem%dealloc(SC, wf%ao%n, wf%n_mo)
      call mem%dealloc(occupation_numbers, wf%n_mo)
!
   end subroutine project_out_cv_and_vc_cuhf
!
!
end module cuhf_class
