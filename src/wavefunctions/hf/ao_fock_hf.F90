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
submodule (hf_class) ao_fock
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
   module subroutine update_fock_and_energy_hf(wf, cumulative)
!!
!!    Wrapper for Update Fock and energy
!!    Written by Tommaso Giovannini, July 2019
!!
!!    Call either cumulative or non-cumulative updating depending on options
!!
!
      use array_utilities, only: copy_and_scale
!
      implicit none
!
      class(hf), intent(inout) :: wf
      logical, intent(in) :: cumulative
      real(dp), dimension(:,:), allocatable :: h
!
      if (cumulative) then
!
         call wf%update_G_cumulative()
!
      else
!
         call wf%update_G_non_cumulative()
!
      endif
!
!
      call copy_and_scale(half, wf%ao_G, wf%ao_fock, wf%ao%n**2)
!
      call mem%alloc(h, wf%ao%n, wf%ao%n)
      call wf%get_ao_h(h)
      call daxpy(wf%ao%n**2, one, h, 1, wf%ao_fock, 1)
!
      if (wf%embedded) then
!
         call wf%embedding%update(wf%ao, wf%ao_density)
         call daxpy(wf%ao%n**2, one, wf%ao%v, 1, wf%ao_fock, 1)
!
      endif
!
!     Energy
!
      wf%energy = wf%calculate_hf_energy_from_fock(wf%ao_fock, h)
!
      call mem%dealloc(h, wf%ao%n, wf%ao%n)
!
   end subroutine update_fock_and_energy_hf
!
!
   module subroutine update_G_cumulative_hf(wf)
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
!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      real(dp), dimension(:, :), allocatable :: G
!
      call output%printf('v', 'Fock matrix construction using density differences')
!
!     Determine density differences
!
      call daxpy(wf%ao%n**2, -one, wf%previous_ao_density, 1, wf%ao_density, 1)
!
      call mem%alloc(G, wf%ao%n, wf%ao%n)
!
      call wf%get_G(wf%ao_density, G)
!
      call daxpy(wf%ao%n**2, one, G, 1, wf%ao_G, 1)
      call mem%dealloc(G, wf%ao%n, wf%ao%n)
!
      call daxpy(wf%ao%n**2, one, wf%previous_ao_density, 1, wf%ao_density, 1)
!
   end subroutine update_G_cumulative_hf
!
!
   module subroutine update_G_non_cumulative_hf(wf)
!!
!!    Update G non cumulative
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Called by solver when a new density has been obtained and
!!    the next Fock and energy is to be computed.
!!
!
!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      call wf%get_G(wf%ao_density, wf%ao_G)
!
   end subroutine update_G_non_cumulative_hf
!
!
   module subroutine get_G_hf(wf, D, G)
!!
!!
!
      use abstract_G_adder_class,      only: abstract_G_adder
      use abstract_G_screener_class,   only: abstract_G_screener
      use ao_G_builder_class,          only: ao_G_builder
!
      use abstract_G_tool_factory_class, only: abstract_G_tool_factory
      use G_tool_factory_class,          only: G_tool_factory
      use J_tool_factory_class,          only: J_tool_factory
      use K_tool_factory_class,          only: K_tool_factory
!
      implicit none
!
      class(hf), intent(inout) :: wf
      real(dp), dimension(wf%ao%n**2, wf%n_densities), intent(in) :: D
      real(dp), dimension(wf%ao%n**2, wf%n_densities), intent(out) :: G
!
      type(timings) :: timer
      real(dp), dimension(:, :), allocatable :: K, J
!
      type(ao_G_builder),             allocatable :: G_builder
      class(abstract_G_screener),     allocatable :: screener
      class(abstract_G_adder),        allocatable :: adder
      class(abstract_G_tool_factory), allocatable :: factory
!
      timer = timings('AO G construction', pl='normal')
      call timer%turn_on()
!
!     Construct the two electron part of the Fock matrix (G),
!     and add the contribution to the Fock matrix
!
      if (wf%coulomb_exchange_separated) then
!
         call mem%alloc(J, wf%ao%n, wf%ao%n)
!
         factory = J_tool_factory(wf%coulomb_threshold)
         call factory%create(screener, adder)
!
         G_builder = ao_G_builder(screener, adder)
         call G_builder%construct(wf%ao, wf%eri_getter, D, J)
!
         call dcopy(wf%ao%n**2, J, 1, G, 1)
!
         call mem%dealloc(J, wf%ao%n, wf%ao%n)
!
         call mem%alloc(K, wf%ao%n, wf%ao%n)
!
         deallocate(screener)
         deallocate(adder)
!
         factory = K_tool_factory(wf%exchange_threshold)
         call factory%create(screener, adder)
!
         G_builder = ao_G_builder(screener, adder)
         call G_builder%construct(wf%ao, wf%eri_getter, D, K)
!
         call daxpy(wf%ao%n**2, one, K, 1, G, 1)
!
         call mem%dealloc(K, wf%ao%n, wf%ao%n)
!
      else
!
         factory = G_tool_factory(wf%coulomb_threshold, wf%exchange_threshold)
         call factory%create(screener, adder)
!
         G_builder = ao_G_builder(screener, adder)
         call G_builder%construct(wf%ao, wf%eri_getter, D, G)
!
      endif
!
      call timer%turn_off()
!
   end subroutine get_G_hf
!
!
   module subroutine construct_ao_G_1der_hf(wf, G_ao, D_ao)
!!
!!    Construct AO G 1der
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    This routine constructs the two-electron part of the Fock matrix G(D),
!!
!!       G_αβqk = sum_γδ 2 * g_αβγδqk D_γδ - sum_γδ g_αδγβqk D_γδ (= G(D)_αβqk),
!!
!!    where contributions from different threads are gathered column blocks
!!    of the incoming F matrix.
!!
!
      use omp_lib
!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%ao%n, wf%ao%n, 3, wf%n_atomic_centers) :: G_ao
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in) :: D_ao
!
      real(dp), dimension((wf%ao%max_sh_size**4)*3*4), target :: g_ABCDqk
      real(dp), dimension(:,:,:,:,:,:), pointer, contiguous :: g_ABCDqk_p
!
      integer :: A, B, C, D, D_max, w, x, y, z, AB, AB_packed
      integer :: w_red, x_red, y_red, z_red, tot_dim, k, q, n_threads, thread
!
      real(dp) :: d1, d2, d3, d4, d5, d6
!
      integer, dimension(4) :: atoms
!
      real(dp) :: deg, deg_AB, deg_CD, deg_AB_CD
!
      real(dp), dimension(3,4) :: temp, temp1, temp2, temp3, temp4, temp5, temp6
!
      real(dp), dimension(:,:,:,:,:), allocatable :: G_ao_t
!
      thread = 0
      n_threads = 1
!$    n_threads = omp_get_max_threads()
      call mem%alloc(G_ao_t, wf%ao%n, wf%ao%n, 3, wf%n_atomic_centers, n_threads)
      G_ao_t = zero
!
!$omp parallel do private(A, B, C, D, D_max, atoms, deg, deg_CD, deg_AB, deg_AB_CD, g_ABCDqk, g_ABCDqk_p, &
!$omp w, x, y, z, w_red, x_red, y_red, z_red, temp, temp1, temp2, temp3, temp4, temp5, temp6, &
!$omp d1, d2, d3, d4, d5, d6, thread, q, k, tot_dim, AB, AB_packed) schedule(dynamic)
      do AB = 1, wf%ao%n_sig_eri_shp
!
!$       thread = omp_get_thread_num()
!
         if (wf%ao%cs_eri_max(AB, 1)*wf%ao%cs_eri_max(1, 1) < wf%coulomb_threshold) cycle
         AB_packed = wf%ao%cs_eri_max_indices(AB, 3)
!
         A = wf%ao%cs_eri_max_indices(AB_packed, 1)
         B = wf%ao%cs_eri_max_indices(AB_packed, 2)
!
         atoms(1) = wf%ao%shell_to_center(A)
         atoms(2) = wf%ao%shell_to_center(B)
!
         deg_AB = real(2-B/A, kind=dp)
!
            do C = 1, A
!
               D_max = (C/A)*B + (1-C/A)*C
               atoms(3) = wf%ao%shell_to_center(C)
!
               do D = 1, D_max
!
                  deg_CD    = real(2-D/C, kind=dp)
                  deg_AB_CD = min(1-C/A+2-min(D/B,B/D), 2)

                  deg = deg_AB*deg_CD*deg_AB_CD ! Shell degeneracy
!
                  atoms(4) = wf%ao%shell_to_center(D)
!
                  call wf%ao%get_eri_1der(g_ABCDqk, A, B, C, D)
!
                  tot_dim = wf%ao%shells(A)%length * &
                            wf%ao%shells(B)%length * &
                            wf%ao%shells(C)%length * &
                            wf%ao%shells(D)%length * 3 * 4
!
                  g_ABCDqk(1:tot_dim) = deg*g_ABCDqk(1:tot_dim)
!
                  g_ABCDqk_p(1 : wf%ao%shells(A)%length, &
                             1 : wf%ao%shells(B)%length, &
                             1 : wf%ao%shells(C)%length, &
                             1 : wf%ao%shells(D)%length, &
                             1 : 3, 1 : 4) => g_ABCDqk(1 : tot_dim)
!
                  do z = wf%ao%shells(D)%first, wf%ao%shells(D)%get_last()
!
                     z_red = z - wf%ao%shells(D)%first + 1
!
                     do y = wf%ao%shells(C)%first, wf%ao%shells(C)%get_last()
!
                        y_red = y - wf%ao%shells(C)%first + 1
!
                        d1 = D_ao(y, z)
!
                        do x = wf%ao%shells(B)%first, wf%ao%shells(B)%get_last()
!
                           x_red = x - wf%ao%shells(B)%first + 1
!
                           d3 = D_ao(x, y)
                           d5 = D_ao(x, z)
!
                           do w = wf%ao%shells(A)%first, wf%ao%shells(A)%get_last()
!
                              d2 = D_ao(w, x)
                              d4 = D_ao(w, y)
                              d6 = D_ao(w, z)
!
                              w_red = w - wf%ao%shells(A)%first + 1
!
                              do k = 1, 4
                                 do q = 1, 3
                                    temp(q, k) = g_ABCDqk_p(w_red, x_red, y_red, z_red, q, k)
                                 enddo
                              enddo
!
                              temp1 = temp*d1
                              temp2 = temp*d2
!
                              temp3 = quarter*temp*d3
                              temp4 = quarter*temp*d4
                              temp5 = quarter*temp*d5
                              temp6 = quarter*temp*d6
!
                              do k = 1, 4
                                 do q = 1, 3
!
                                    G_ao_t(w, x, q, atoms(k), thread+1) = &
                                          G_ao_t(w, x, q, atoms(k), thread+1) + temp1(q, k)
                                    G_ao_t(y, x, q, atoms(k), thread+1) = &
                                          G_ao_t(y, x, q, atoms(k), thread+1) - temp6(q, k)
                                    G_ao_t(y, z, q, atoms(k), thread+1) = &
                                          G_ao_t(y, z, q, atoms(k), thread+1) + temp2(q, k)
                                    G_ao_t(w, z, q, atoms(k), thread+1) = &
                                          G_ao_t(w, z, q, atoms(k), thread+1) - temp3(q, k)
                                    G_ao_t(x, z, q, atoms(k), thread+1) = &
                                          G_ao_t(x, z, q, atoms(k), thread+1) - temp4(q, k)
                                    G_ao_t(w, y, q, atoms(k), thread+1) = &
                                          G_ao_t(w, y, q, atoms(k), thread+1) - temp5(q, k)
!
                                 enddo
                              enddo
!
                           enddo
                        enddo
                     enddo
                  enddo
!
               enddo
            enddo
         enddo
!$omp end parallel do
!
      G_ao = zero
!
      do thread = 1, n_threads
!
         call daxpy(3*wf%n_atomic_centers*wf%ao%n**2, one, G_ao_t(1,1,1,1,thread), 1, G_ao, 1)
!
      enddo
!
      call mem%dealloc(G_ao_t, wf%ao%n, wf%ao%n, 3, wf%n_atomic_centers, n_threads)
!
   end subroutine construct_ao_G_1der_hf
!
!
   module subroutine set_screening_and_precision_thresholds_hf(wf, gradient_threshold)
!!
!!    Set screening and precision thresholds
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the screening thresholds for Coulomb and exchange
!!    integrals given the convergence threshold for the gradient
!!    unless other thresholds are already set on input.
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      real(dp), intent(in) :: gradient_threshold
!
      real(dp) :: epsilon_ ! Libint ERI precision
!
      real(dp) :: integral_cutoff
!
!     If not specified by user, set Coulomb and exchange thresholds
!     (For ref., this will give 10-12 and 10-10 for gradient threshold: 10-6)
!
      if (input%is_keyword_present('coulomb threshold', 'solver scf')) then
!
         call input%get_keyword('coulomb threshold',  &
                                        'solver scf', &
                                        wf%coulomb_threshold)
!
      else
!
         wf%coulomb_threshold  = gradient_threshold*1.0d-6
!
      endif
!
      if (input%is_keyword_present('exchange threshold', 'solver scf')) then
!
         call input%get_keyword('exchange threshold', &
                                        'solver scf', &
                                        wf%exchange_threshold)
!
      else
!
         wf%exchange_threshold = gradient_threshold*1.0d-4
!
      endif
!
!     Electron repulsion integrals
!
      if (input%is_keyword_present('integral cutoff', 'solver scf')) then
!
         call input%get_keyword('integral cutoff',          &
                                          'solver scf',     &
                                          integral_cutoff)
!
      else
!
         integral_cutoff = wf%coulomb_threshold
!
      endif
!
      call wf%ao%set_eri_cutoff(integral_cutoff)
!
!     If not requested by user, set Libint integral accuracy according to thresholds
!
      if (input%is_keyword_present('integral precision', 'solver scf')) then
!
         call input%get_keyword('integral precision',       &
                                          'solver scf',     &
                                          epsilon_)
!
!        Warn user about tampering with integral precision (Libint epsilon)
!
         call output%warning_msg("Setting a specific integral precision can cause instabilities.&
                                 & Be aware that the 'precision' given to Libint is not exact&
                                 & but approximate. Don't use unless you know what you are doing.")
!
!
      else
!
         epsilon_ = integral_cutoff**2
!
      endif
!
      call wf%ao%set_libint_epsilon(epsilon_)
!
!     One electron integrals
!
      if (input%is_keyword_present('one-electron integral cutoff', 'solver scf')) then
!
         call input%get_keyword('one-electron integral cutoff',   &
                                           'solver scf',          &
                                           integral_cutoff)
!
         call wf%ao%set_oei_cutoff(integral_cutoff)
!
      else
!
         call wf%ao%set_oei_cutoff(wf%coulomb_threshold*1.0d-5)
!
      endif
!
!     Currently, exchange threshold must be >= coulomb threshold in the Fock construction
!
      if (wf%exchange_threshold .lt. wf%coulomb_threshold) then
!
         call output%printf('m', 'Exchange threshold is restricted to being higher than, or&
                                 & equal, to the Coulomb threshold. Setting them equal to&
                                 & each other.')
!
         wf%exchange_threshold = wf%coulomb_threshold
!
      endif
!
   end subroutine set_screening_and_precision_thresholds_hf
!
!
   module subroutine set_coulomb_exchange_separation_hf(wf)
!!
!!    Set Coulomb and exchange separation
!!    Written by Linda Goletto, Aug 2020
!!
!!    If the 'coulomb exchange terms' keyword is defined in the input,
!!    sets the first element of the wf%coulomb_exchange_separation list
!!    to 'requested' and the second element to the input request;
!!    if the 'coulomb exchange terms' keyword is not defined in the input,
!!    sets the first element of the wf%coulomb_exchange_separation list
!!    to 'default' and the second element to either 'separated' or
!!    'collective', according to the number of atoms being larger or
!!    smaller/equal to a limit set to 200.
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      integer :: n_atoms_limit = 200
      character(len=200) :: coulomb_exchange_separation
!
      if (input%is_keyword_present('coulomb exchange terms', 'solver scf')) then
!
         call input%get_keyword('coulomb exchange terms', &
                                           'solver scf',  &
                                            coulomb_exchange_separation)
!
         wf%coulomb_exchange_separated = (coulomb_exchange_separation == 'separated')
!
      else
!
         if (wf%n_atomic_centers .le. n_atoms_limit) then
!
            wf%coulomb_exchange_separated = .false.
!
         else
!
            wf%coulomb_exchange_separated = .true.
!
         endif
!
      endif
!
      if (wf%coulomb_exchange_separated) then
!
         call output%printf('v', 'Will perform Coulomb and exchange terms in the Fock&
                                                   & matrix separately', fs='(/t3,a)')
!
      else
!
         call output%printf('v', 'Will perform Coulomb and exchange terms in the Fock&
                                                 & matrix collectively', fs='(/t3,a)')
!
      endif
!
   end subroutine set_coulomb_exchange_separation_hf
!
!
   module subroutine set_gradient_threshold_hf(wf, gradient_threshold)
!!
!!    Set gradient thresholds
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the screening thresholds for Coulomb and exchange
!!    integrals given the convergence threshold for the gradient
!!    unless other thresholds are already set on input.
      implicit none
!
      class(hf), intent(inout) :: wf
      real(dp), intent(in) :: gradient_threshold
!
      wf%gradient_threshold = gradient_threshold
!
      call wf%set_screening_and_precision_thresholds(wf%gradient_threshold)
!
   end subroutine set_gradient_threshold_hf
!
!
end submodule ao_fock

