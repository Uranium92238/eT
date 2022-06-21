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
submodule (uhf_class) ao_fock
!
!!
!!    AO Fock submodule
!!
!!    Collects the routines used in the construction of the AO Fock matrix.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine update_fock_and_energy_cumulative_uhf(wf, prev_ao_density)
!!
!!    Update Fock and energy cumulatively
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    This routine guides the construction of the Fock matrix (or matrices for
!!    unrestricted wavefunctions) from the current AO density (or densities).
!!    It is called by the solver and is overwritten for unrestricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
      real(dp), dimension(wf%ao%n**2, wf%n_densities), intent(in) :: prev_ao_density
      real(dp), dimension(:,:), allocatable :: h
!
      logical :: cumulative
!
      call output%printf('v', 'Fock matrix construction using density differences')
!
!     Here, the previous AO density is sent as [D_a D_b],
!     where each is full square
!
      call daxpy(wf%ao%n**2, -one, prev_ao_density, 1, wf%ao_density_a, 1)
      call daxpy(wf%ao%n**2, -one, prev_ao_density(1, 2), 1, wf%ao_density_b, 1)
!
      call daxpy(wf%ao%n**2, -one, prev_ao_density, 1, wf%ao_density, 1)
      call daxpy(wf%ao%n**2, -one, prev_ao_density(1, 2), 1, wf%ao_density, 1)
!
      cumulative = .true.
!
      call mem%alloc(h, wf%ao%n, wf%ao%n)
      call wf%get_ao_h(h)
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_a, wf%ao_density_b, h, cumulative)
      call mem%dealloc(h, wf%ao%n, wf%ao%n)
!
!     Restore density
!
      call daxpy(wf%ao%n**2, one, prev_ao_density, 1, wf%ao_density_a, 1)
      call daxpy(wf%ao%n**2, one, prev_ao_density(1, 2), 1, wf%ao_density_b, 1)
!
      call daxpy(wf%ao%n**2, one, prev_ao_density, 1, wf%ao_density, 1)
      call daxpy(wf%ao%n**2, one, prev_ao_density(1, 2), 1, wf%ao_density, 1)
!
   end subroutine update_fock_and_energy_cumulative_uhf
!
!
   module subroutine update_fock_and_energy_non_cumulative_uhf(wf)
!!
!!    Update Fock and energy
!!    Written by Eirik F. Kjønstad, Sep 2018
!!    Modified for QM/MM by Tommaso Giovannini, July 2019
!!
!!    This routine guides the construction of the Fock matrix (or matrices for
!!    unrestricted wavefunctions) from the current AO density (or densities).
!!    It is called by the solver and is overwritten for unrestricted
!!    wavefunctions.
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
      real(dp), dimension(:,:), allocatable :: h
      real(dp), dimension(:,:), allocatable :: h_wx_eff
!
!     QM/MM and PCM specific work
!
      call mem%alloc(h, wf%ao%n, wf%ao%n)
      call wf%get_ao_h(h)
!
      if (wf%embedded) then
!
         call wf%embedding%update(wf%ao, wf%ao_density)
!
         call mem%alloc(h_wx_eff, wf%ao%n, wf%ao%n)
!
         call dcopy(wf%ao%n**2, h, 1, h_wx_eff, 1)
         call daxpy(wf%ao%n**2, one, wf%ao%v, 1, h_wx_eff, 1)
!
         call wf%construct_ao_spin_fock(wf%ao_density,   &
                                       wf%ao_density_a,  &
                                       wf%ao_density_b,  &
                                       h_wx_eff, cumulative=.false.)
!
         call mem%dealloc(h_wx_eff, wf%ao%n, wf%ao%n)
!
      else
!
         call wf%construct_ao_spin_fock(wf%ao_density,   &
                                        wf%ao_density_a, &
                                        wf%ao_density_b, &
                                        h, cumulative=.false.)
!
      endif
!
      call mem%dealloc(h, wf%ao%n, wf%ao%n)
!
   end subroutine update_fock_and_energy_non_cumulative_uhf
!
!
   module subroutine update_fock_and_energy_uhf(wf, cumulative)
!!
!!    Update Fock and energy wrapper
!!    Written by Tommaso Giovannini, July 2019
!!
!!    Wrapper for cumulative or non-cumulative subroutines
!!    depending on the path
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
      logical, intent(in) :: cumulative
      real(dp), dimension(:,:), allocatable :: h
!
      if (.not. cumulative) then
!
         call wf%update_fock_and_energy_non_cumulative()
!
      else
!
         if (wf%embedded) then
!
            call wf%update_fock_and_energy_non_cumulative()
!
         else
!
            call wf%update_fock_and_energy_cumulative(wf%previous_ao_density)
!
         endif
!
      endif
!
      call mem%alloc(h, wf%ao%n, wf%ao%n)
      call wf%get_ao_h(h)
      call wf%calculate_uhf_energy(h)
      call mem%dealloc(h, wf%ao%n, wf%ao%n)
!      
   end subroutine update_fock_and_energy_uhf
!
!
   module subroutine construct_ao_spin_fock_uhf(wf, D, D_alpha, D_beta, h_wx, cumulative)
!!
!!    Construct AO spin Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    The routine computes the alpha and beta Fock matrix:
!!
!!       F_αβ^alpha = h_αβ + sum_γδ g_αβγδ D_γδ - sum_γδ g_αδγβ D_γδ^alpha
!!       F_αβ^beta  = h_αβ + sum_γδ g_αβγδ D_γδ - sum_γδ g_αδγβ D_γδ^beta
!!
!!    Here the superscript refers to the spin function, while the subscripts
!!    are AO indices. In contrast to the restricted routine, this one does
!!    not calculate the energy - a separate call is required to get the
!!    unrestricted Hartree-Fock energy.
!!
!
!
      use abstract_G_adder_class,      only: abstract_G_adder
      use abstract_G_screener_class,   only: abstract_G_screener
      use ao_G_builder_class,          only: ao_G_builder
!
      use abstract_G_tool_factory_class, only: abstract_G_tool_factory
      use J_tool_factory_class,          only: J_tool_factory
      use K_tool_factory_class,          only: K_tool_factory
!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D_alpha, D_beta
!
      class(ao_G_builder),              allocatable :: G_builder
      class(abstract_G_screener),      allocatable :: screener
      class(abstract_G_adder),         allocatable :: adder
      class(abstract_G_tool_factory),  allocatable :: factory
!
      logical, intent(in), optional :: cumulative
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: h_wx
!
      logical :: local_cumulative
!
      real(dp), dimension(:,:), allocatable :: X_alpha, X_beta, J
!
      type(timings), allocatable :: timer
!
      timer = timings('AO Fock construction', 'normal')
      call timer%turn_on()
!
!     Set thresholds to ignore Coulomb and exchange terms,
!     as well as the desired Libint integral precision,
!     and determine whether the construction should be
!     cumulative or not
!
      local_cumulative = .false.
      if (present(cumulative)) local_cumulative = cumulative
!
      call mem%alloc(J, wf%ao%n, wf%ao%n)
!
      factory = J_tool_factory(wf%coulomb_threshold)
      call factory%create(screener, adder)
!
      G_builder = ao_G_builder(screener, adder)
      call G_builder%construct(wf%ao, wf%eri_getter, D, J)
!
      call mem%alloc(X_alpha, wf%ao%n, wf%ao%n)
!
      deallocate(screener)
      deallocate(adder)
!
      factory = K_tool_factory(wf%exchange_threshold)
      call factory%create(screener, adder)
      deallocate(factory)
!
      G_builder = ao_G_builder(screener, adder)
      call G_builder%construct(wf%ao, wf%eri_getter, D_alpha, X_alpha)
      deallocate(G_builder)
!
      call mem%alloc(X_beta, wf%ao%n, wf%ao%n)
!
      G_builder = ao_G_builder(screener, adder)
      call G_builder%construct(wf%ao, wf%eri_getter, D_beta, X_beta)
      deallocate(G_builder)
!
!     Add the accumulated Fock matrix F into the correct Fock matrix
!     (i.e., either the alpha or beta Fock matrix )
!
      if (.not. local_cumulative) call dcopy(wf%ao%n**2, h_wx, 1, wf%ao_fock_a, 1)
!
      call daxpy(wf%ao%n**2, half, J, 1, wf%ao_fock_a, 1)
      call daxpy(wf%ao%n**2, one, X_alpha, 1, wf%ao_fock_a, 1)
!
      if (.not. local_cumulative) call dcopy(wf%ao%n**2, h_wx, 1, wf%ao_fock_b, 1)
!
      call daxpy(wf%ao%n**2, half, J, 1, wf%ao_fock_b, 1)
      call daxpy(wf%ao%n**2, one, X_beta, 1, wf%ao_fock_b, 1)
!
      call mem%dealloc(J, wf%ao%n, wf%ao%n)
      call mem%dealloc(X_alpha, wf%ao%n, wf%ao%n)
      call mem%dealloc(X_beta, wf%ao%n, wf%ao%n)
      call timer%turn_off()
!
   end subroutine construct_ao_spin_fock_uhf
!
end submodule ao_fock
