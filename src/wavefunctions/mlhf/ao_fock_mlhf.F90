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
submodule (mlhf_class) ao_fock
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
   module subroutine update_G_cumulative_mlhf(wf)
!!
!!    Update G cumulative
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Called by solver when a new density has been obtained and
!!    the next Fock and energy is to be computed. By cumulatively
!!    we mean using the density change to build the Fock matrix
!!    in the iterative loop.
!!
!!    Modified by SDF, Sep 2020:
!!
!!    Cumulative construction of two-electron part and not
!!    full Fock matrix
!!
      implicit none
!
      class(mlhf), intent(inout) :: wf
!
      real(dp), dimension(:,:), allocatable :: G
!
      call output%printf('v', 'Fock matrix construction using density differences')
!
      call daxpy(wf%ao%n**2, -one, wf%previous_ao_density, 1, wf%ao_density, 1)
!
      call mem%alloc(G, wf%ao%n, wf%ao%n)
!
      if (wf%C_screening) then
!
         call wf%get_G_MO_screened(wf%ao_density, G)
!
      else
!
         call wf%get_G(wf%ao_density, G)
!
      endif
!
      call daxpy(wf%ao%n**2, half, G, 1, wf%ao_fock, 1)
      call mem%dealloc(G, wf%ao%n, wf%ao%n)
!
      call daxpy(wf%ao%n**2, one, wf%previous_ao_density, 1, wf%ao_density, 1)
!
   end subroutine update_G_cumulative_mlhf
!
!
   module subroutine update_G_non_cumulative_mlhf(wf)
!!
!!    Update G non cumulative
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Called by solver when a new density has been obtained and
!!    the next Fock and energy is to be computed.
!!
!
      implicit none
!
      class(mlhf), intent(inout) :: wf
!
      if (wf%C_screening) then
!
         call wf%get_G_MO_screened(wf%ao_density, wf%ao_fock)
!
      else
!
         call wf%get_G(wf%ao_density, wf%ao_fock)
!
      endif
!
      call dscal(wf%ao%n**2, half, wf%ao_fock, 1)
!
   end subroutine update_G_non_cumulative_mlhf
!
!
   module subroutine get_G_MO_screened_mlhf(wf, D, G)
!!
!!
!
      use abstract_G_adder_class,         only: abstract_G_adder
      use abstract_G_tool_factory_class,  only: abstract_G_tool_factory
      use abstract_G_screener_class,      only: abstract_G_screener
      use ao_G_builder_class,             only: ao_G_builder
!
      use G_tool_MO_screening_factory_class, only: G_tool_MO_screening_factory
      use J_tool_MO_screening_factory_class, only: J_tool_MO_screening_factory
      use K_tool_MO_screening_factory_class, only: K_tool_MO_screening_factory
!
      implicit none
!
      class(mlhf), intent(inout) :: wf
!
      real(dp), dimension(wf%ao%n**2, wf%n_densities), intent(in) :: D
      real(dp), dimension(wf%ao%n**2, wf%n_densities), intent(out) :: G
!
      type(timings) :: timer
      real(dp), dimension(:, :), allocatable :: K, J
!
      type(ao_G_builder),              allocatable :: G_builder
      class(abstract_G_screener),      allocatable :: screener
      class(abstract_G_adder),         allocatable :: adder
      class(abstract_G_tool_factory),  allocatable :: factory
!
      timer = timings('AO G construction', pl='normal')
      call timer%turn_on()
!
      if (wf%coulomb_exchange_separated) then
!
         call mem%alloc(J, wf%ao%n, wf%ao%n)
!
         factory = J_tool_MO_screening_factory(wf%coulomb_threshold, &
                                               wf%ao%n, wf%n_mo, wf%orbital_coefficients)
         call factory%create(screener, adder)
!
         G_builder = ao_G_builder(screener, adder)
         call G_builder%construct(wf%ao, wf%eri_getter, D, J)
!
         call dcopy(wf%ao%n**2, J, 1, G, 1)
!
         call mem%dealloc(J, wf%ao%n, wf%ao%n)
         call mem%alloc(K, wf%ao%n, wf%ao%n)
!
         deallocate(screener)
         deallocate(adder)
!
         factory = K_tool_MO_screening_factory(wf%exchange_threshold, &
                                               wf%ao%n, wf%n_mo, wf%orbital_coefficients)
         call factory%create(screener, adder)
!
         G_builder = ao_G_builder(screener, adder)
         call G_builder%construct(wf%ao, wf%eri_getter, D, K)
!
         call daxpy(wf%ao%n**2, one, K, 1, G, 1)
         call mem%dealloc(K, wf%ao%n, wf%ao%n)
!
      else
!
         factory = G_tool_MO_screening_factory(wf%coulomb_threshold, wf%exchange_threshold, &
                                               wf%ao%n, wf%n_mo, wf%orbital_coefficients)
         call factory%create(screener, adder)
!
         G_builder = ao_G_builder(screener, adder)
         call G_builder%construct(wf%ao, wf%eri_getter, D, G)
!
      endif
!
      call timer%turn_off()
!
   end subroutine get_G_MO_screened_mlhf
!
!
end submodule ao_fock

