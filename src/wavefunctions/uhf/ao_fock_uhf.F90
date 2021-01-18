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
   module subroutine construct_mo_fock_uhf(wf)
!!
!!    Construct MO Fock
!!    Written by Eirik F. Kjønstad, Mar 2019
!!
!!    Give notice to user that it does not yet exist.
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
      call output%printf('m', 'Requested MO transformation of Fock matrix, but &
                         &this is not yet implemented for (a0).', &
                         chars=[trim(wf%name_)], ffs='(/t3,a)')
!
   end subroutine construct_mo_fock_uhf
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
      class(uhf) :: wf
!
      real(dp), dimension(wf%n_ao**2, wf%n_densities), intent(in) :: prev_ao_density
!
      logical :: cumulative
!
!     Here, the previous AO density is sent as [D_a D_b],
!     where each is full square
!
      call daxpy(wf%n_ao**2, -one, prev_ao_density, 1, wf%ao_density_a, 1)
      call daxpy(wf%n_ao**2, -one, prev_ao_density(1, 2), 1, wf%ao_density_b, 1)
!
      call daxpy(wf%n_ao**2, -one, prev_ao_density, 1, wf%ao_density, 1)
      call daxpy(wf%n_ao**2, -one, prev_ao_density(1, 2), 1, wf%ao_density, 1)
!
      cumulative = .true.
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_a, 'alpha', &
                      wf%ao_h, cumulative)
!
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_b, 'beta', &
                      wf%ao_h, cumulative)
!
      call daxpy(wf%n_ao**2, one, prev_ao_density, 1, wf%ao_density_a, 1)
      call daxpy(wf%n_ao**2, one, prev_ao_density(1, 2), 1, wf%ao_density_b, 1)
!
      call daxpy(wf%n_ao**2, one, prev_ao_density, 1, wf%ao_density, 1)
      call daxpy(wf%n_ao**2, one, prev_ao_density(1, 2), 1, wf%ao_density, 1)
!
      call wf%calculate_uhf_energy(wf%ao_h)
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
      class(uhf) :: wf
!      
      real(dp), dimension(:, :), allocatable            :: h_wx_eff
!
      if(wf%system%mm_calculation .and. wf%system%mm%forcefield .eq. 'non-polarizable') then

         call mem%alloc(h_wx_eff, wf%n_ao, wf%n_ao)
!
         call dcopy(wf%n_ao**2, wf%ao_h, 1, h_wx_eff, 1)
         call daxpy(wf%n_ao**2, half, wf%nopol_h_wx, 1, h_wx_eff, 1)
!
         call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_a, 'alpha', h_wx_eff)
         call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_b, 'beta', h_wx_eff)
!
         call mem%dealloc(h_wx_eff, wf%n_ao, wf%n_ao) 
!
      else

         call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_a, 'alpha', wf%ao_h)
         call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_b, 'beta', wf%ao_h)
!
      endif
!      
      if(wf%system%mm_calculation .and. wf%system%mm%forcefield .ne. 'non-polarizable') &
         call wf%update_fock_mm()
!         
      if(wf%system%pcm_calculation) call wf%update_fock_pcm()
!      
      call wf%calculate_uhf_energy(wf%ao_h)
!
   end subroutine update_fock_and_energy_non_cumulative_uhf
!
!
   module subroutine update_fock_and_energy_uhf(wf, prev_ao_density)
!!
!!    Update Fock and energy wrapper
!!    Written by Tommaso Giovannini, July 2019
!!
!!    Wrapper for cumulative or non-cumulative subroutines
!!    depending on the path
!!
      implicit none
!
      class(uhf) :: wf
!
      real(dp), dimension(wf%n_ao**2, wf%n_densities), intent(in), optional :: prev_ao_density
!
      if (.not. present(prev_ao_density)) then 
!
          call wf%update_fock_and_energy_non_cumulative()
!
      else 
!      
         if (.not. wf%system%mm_calculation .and. .not. wf%system%pcm_calculation) then
!         
            call wf%update_fock_and_energy_cumulative(prev_ao_density)
!            
         else
!         
            if (wf%system%mm%forcefield .eq. 'non-polarizable') then
!
               call wf%update_fock_and_energy_cumulative(prev_ao_density)
!
            else 
!         
               call wf%update_fock_and_energy_non_cumulative()
!               
            endif
!         
         endif
!
      endif
!      
   end subroutine update_fock_and_energy_uhf
!
!
   module subroutine construct_ao_spin_fock_uhf(wf, D, D_sigma, sigma, &
                     h_wx, cumulative)
!!
!!    Construct AO spin Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    The routine computes the alpha or beta Fock matrix, depending
!!    on the value of the spin 'sigma' (='alpha' or 'beta'):
!!
!!       F_αβ^alpha = h_αβ + sum_γδ g_αβγδ D_γδ - sum_γδ g_αδγβ D_γδ^alpha
!!       F_αβ^beta  = h_αβ + sum_γδ g_αβγδ D_γδ - sum_γδ g_αδγβ D_γδ^beta
!!
!!    Here the superscript refers to the spin function, while the subscripts
!!    are AO indices. In contrast to the restricted routine, this one does
!!    not calculate the energy - a separate call is required to get the
!!    unrestricted Hartree-Fock energy.
!!
      implicit none
!
      class(uhf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D_sigma
!
      character(len=*), intent(in) :: sigma
!
      logical, intent(in), optional :: cumulative
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      integer :: thread = 0, n_threads = 1
      logical :: local_cumulative
!
      real(dp), dimension(:,:), allocatable :: F, shp_density_schwarz
!
      real(dp) :: coulomb_thr, exchange_thr, precision_thr    ! Actual thresholds
!
      integer :: n_sig_shp
!
      real(dp) :: max_D_schwarz, max_eri_schwarz
!
      real(dp), dimension(:,:), allocatable :: scaled_D_sigma ! = 2 * D_sigma
!
      type(timings), allocatable :: timer 
!
      timer = timings('AO Fock construction (' // trim(sigma) // ')', 'normal')
      call timer%turn_on()
!
!     Set thresholds to ignore Coulomb and exchange terms,
!     as well as the desired Libint integral precision,
!     and determine whether the construction should be
!     cumulative or not
!
      coulomb_thr   = wf%coulomb_threshold
      exchange_thr  = wf%exchange_threshold
      precision_thr = wf%libint_epsilon
!
      local_cumulative = .false.
      if (present(cumulative)) then
!
         if (cumulative) then
!
            local_cumulative = .true.
!
         else
!
            local_cumulative = .false.
!
         endif
!
      endif
!
!     Compute number of significant ERI shell pairs (the Fock construction
!     only loops over these shell pairs) and the maximum element
!
      call wf%get_n_sig_eri_shp(n_sig_shp)
      max_eri_schwarz = get_abs_max(wf%shp_eri_schwarz, wf%system%n_s*(wf%system%n_s + 1)/2)
!
!     Construct the Coulomb two electron part of the Fock matrix, using the screening vectors
!     and parallellizing over available threads (each gets its own copy of the Fock matrix)
!
      call mem%alloc(shp_density_schwarz, wf%system%n_s, wf%system%n_s)
      call wf%construct_shp_density_schwarz(shp_density_schwarz, D)
      max_D_schwarz = get_abs_max(shp_density_schwarz, wf%system%n_s**2)
!
!$      n_threads = omp_get_max_threads()
!
      call mem%alloc(F, wf%n_ao, wf%n_ao*n_threads) ! [F(thread 1) F(thread 2) ...]
      F = zero
!
      call wf%construct_coulomb_ao_G(F, D, n_threads, max_D_schwarz, max_eri_schwarz,         &
                                                shp_density_schwarz,                          &
                                                n_sig_shp, coulomb_thr, precision_thr,         &
                                                wf%system%shell_limits)
!
!     Construct the Coulomb two electron part of the Fock matrix, using the screening vectors
!     and parallellizing over available threads (each gets its own copy of the Fock matrix)
!
      call mem%alloc(scaled_D_sigma, wf%n_ao, wf%n_ao)
      scaled_D_sigma = two*D_sigma
!
      call wf%construct_shp_density_schwarz(shp_density_schwarz, scaled_D_sigma)
      max_D_schwarz = get_abs_max(shp_density_schwarz, wf%system%n_s**2)
!
      call wf%construct_exchange_ao_G(F, scaled_D_sigma, n_threads, max_D_schwarz,    &
                                      max_eri_schwarz, shp_density_schwarz, n_sig_shp, &
                                      exchange_thr, precision_thr, wf%system%shell_limits)
!
      call mem%dealloc(shp_density_schwarz, wf%system%n_s, wf%system%n_s)
      call mem%dealloc(scaled_D_sigma, wf%n_ao, wf%n_ao)
!
!     Add the accumulated Fock matrix F into the correct Fock matrix
!     (i.e., either the alpha or beta Fock matrix )
!
      if (trim(sigma) == 'alpha') then
!
         if (.not. local_cumulative) wf%ao_fock_a = zero
         do thread = 1, n_threads
!
            call daxpy(wf%n_ao**2, one, F(1, (thread-1)*wf%n_ao + 1), 1, wf%ao_fock_a, 1)
!
         enddo
!
         call symmetric_sum(wf%ao_fock_a, wf%n_ao)
         wf%ao_fock_a = wf%ao_fock_a*half
!
         if (.not. local_cumulative) wf%ao_fock_a = wf%ao_fock_a + h_wx
!
      elseif (trim(sigma) == 'beta') then
!
         if (.not. local_cumulative) wf%ao_fock_b = zero
         do thread = 1, n_threads
!
            call daxpy(wf%n_ao**2, one, F(1, (thread-1)*wf%n_ao + 1), 1, wf%ao_fock_b, 1)
!
         enddo
!
         call symmetric_sum(wf%ao_fock_b, wf%n_ao)
         wf%ao_fock_b = wf%ao_fock_b*half
!
         if (.not. local_cumulative) wf%ao_fock_b = wf%ao_fock_b + h_wx
!
      else
!
         call output%error_msg('Did not recognize spin variable in construct_ao_fock:' // trim(sigma))
!
      endif
!
      call mem%dealloc(F, wf%n_ao, wf%n_ao*n_threads)
      call timer%turn_off()
!
   end subroutine construct_ao_spin_fock_uhf
!
!
   module subroutine update_fock_mm_uhf(wf)
!!
!!    Update alpha and beta Fock matrices with polarizable QM/MM terms
!!    For now: QM/FQ model (see mm_class and output file)
!!    Written by Tommaso Giovannini, July 2019 for QM/MM
!!
      implicit none
!
      class(uhf) :: wf
!
      real(dp), dimension(:),    allocatable                  :: potential_points
      integer :: i
!     
      if(wf%system%mm%forcefield.eq.'fq') then
!      
         if(.not.allocated(potential_points)) call mem%alloc(potential_points, wf%system%mm%n_atoms)
!
         call zero_array(wf%pol_emb_fock,wf%n_ao*wf%n_ao)
         call zero_array(wf%system%mm%pol_emb_rhs,wf%system%mm%n_variables)
!         
!        electrostatic potential contracted with density : \sum_i V_mu(D_mu)(r_i)
!
         call wf%construct_ao_electrostatics(0,1,'prop',wf%system%mm%n_atoms,wf%system%mm%coordinates, &
                                             property_points=potential_points,ao_density=wf%ao_density)
!
!        rhs for fq: -chi - V(D)
!
         wf%system%mm%pol_emb_rhs(1:wf%system%mm%n_atoms) = -wf%system%mm%chi + potential_points
!         
!        solve q=D^-1 (-chi - V(D)) 
! 
         call dgemm('N', 'N',                   &
                     wf%system%mm%n_variables,  &
                     1,                         &
                     wf%system%mm%n_variables,  &
                     one,                       &
                     wf%system%mm%fq_matrix,    &
                     wf%system%mm%n_variables,  &
                     wf%system%mm%pol_emb_rhs,         &
                     wf%system%mm%n_variables,  &
                     zero,                      &
                     wf%system%mm%pol_emb_lhs,            &
                     wf%system%mm%n_variables)
!
!
         call output%print_separator('verbose', 67, fs='(/t3,a)')
!
         call output%printf('v', 'Atom          FQ LHS             FQ RHS       &
                            & QM Potential@FQs', fs='(t6,a)')
!            
         do i = 1, wf%system%mm%n_atoms
!           
            call output%printf('v', '(i4)      (e13.6)      (e13.6)      (e13.6)', &
                               ints=[i], reals=[wf%system%mm%pol_emb_lhs(i), &
                               wf%system%mm%pol_emb_rhs(i), &
                               potential_points(i)], fs='(t6,a)')
!           
         enddo
!           
         call output%print_separator('verbose', 67)
!
!
!        put FQ charges into charge (I am discrading langrangian multipliers)
!
         wf%system%mm%charge = wf%system%mm%pol_emb_lhs(1:wf%system%mm%n_atoms)
!
!        Fock creation: F_munu = \sum_i q_i V_munu(r_i)
!
         call wf%construct_ao_electrostatics(0,0,'fock',wf%system%mm%n_atoms,wf%system%mm%coordinates, &
                                             elec_fock=wf%pol_emb_fock,charges=wf%system%mm%charge) 
!
         wf%ao_fock_a = wf%ao_fock_a + half*wf%pol_emb_fock
         wf%ao_fock_b = wf%ao_fock_b + half*wf%pol_emb_fock
!
!
         call output%print_matrix('debug', 'Total QM Density', wf%ao_density, wf%n_ao, wf%n_ao)
!
         call output%print_matrix('debug', 'QM Density - Spin alpha', &
                                  wf%ao_density_a, wf%n_ao, wf%n_ao)
!
         call output%print_matrix('debug', 'QM Density - Spin beta', &
                                  wf%ao_density_b, wf%n_ao, wf%n_ao)
!
         call output%print_matrix('debug', 'FQ Fock - Spin alpha + beta', &
                                  wf%pol_emb_fock, wf%n_ao, wf%n_ao)
!
         call output%print_matrix('debug', 'QM/FQ Fock - Spin alpha', wf%ao_fock_a, wf%n_ao, wf%n_ao)
!
         call output%print_matrix('debug', 'QM/FQ Fock - Spin beta', wf%ao_fock_b, wf%n_ao, wf%n_ao)
!
!
         call mem%dealloc(potential_points, wf%system%mm%n_atoms)
!         
      else
!      
         call output%error_msg('The only available polarizable force field is fq')
!         
      endif
!
!
   end subroutine update_fock_mm_uhf
!
!
   module subroutine update_fock_pcm_uhf(wf)
!!
!!    Update Fock PCM 
!!    Written by Tommaso Giovannini, Oct 2019 
!!
!!    The QM Fock is updated with the contributions coming 
!!    from the PCM:
!!       q*V_wx
!!
!!    Done by interfacing to PCMSolver
!!
      implicit none
!
      class(uhf) :: wf
!
      integer :: i
! 
      call zero_array(wf%pcm_fock,wf%n_ao*wf%n_ao)
      call zero_array(wf%system%pcm%pcm_rhs,wf%system%pcm%n_tesserae)
!      
!     electrostatic potential contracted with density : \sum_i V_mu(D_mu)(r_i)
!
      call wf%construct_ao_electrostatics(0, 1, 'prop', &
                                          wf%system%pcm%n_tesserae, &
                                          wf%system%pcm%grid_coord*bohr_to_angstrom, &
                                          property_points=wf%system%pcm%pcm_rhs, &
                                          ao_density=wf%ao_density) 
!      
!     solve q=D^-1 (V(D)) 
! 
      call wf%system%pcm%set_surface_function('NucMEP')
!                                          
      call wf%system%pcm%compute_asc('NucMEP', 'NucASC')
!                                          
      call wf%system%pcm%get_surface_function('NucASC')
!
      call output%print_separator('verbose', 67, fs='(/t3,a)')
!
      call output%printf('v', 'Atom         PCM ASC            PCM RHS', fs='(t6,a)')
!
      do i = 1, wf%system%pcm%n_tesserae
!
         call output%printf('v', '(i4)      (e13.6)      (e13.6)', ints=[i], &
                            reals=[wf%system%pcm%charges(i), &
                            wf%system%pcm%pcm_rhs(i)], fs='(t6,a)')
!
      enddo
!
      call output%print_separator('verbose', 67)
!
!     Fock creation: F_munu = \sum_i q_i V_munu(r_i)
!
      call wf%construct_ao_electrostatics(0, 0, 'fock', &
                                          wf%system%pcm%n_tesserae, &
                                          wf%system%pcm%grid_coord*bohr_to_angstrom, & 
                                          elec_fock=wf%pcm_fock, &
                                          charges=wf%system%pcm%charges) 
!
      wf%ao_fock_a = wf%ao_fock_a + half*wf%pcm_fock
      wf%ao_fock_b = wf%ao_fock_b + half*wf%pcm_fock
!
      call output%print_matrix('debug', 'Total QM Density', wf%ao_density, wf%n_ao, wf%n_ao)
!
      call output%print_matrix('debug', 'QM Density - Spin alpha', &
                               wf%ao_density_a, wf%n_ao, wf%n_ao)
!
      call output%print_matrix('debug', 'QM Density - Spin beta', &
                               wf%ao_density_b, wf%n_ao, wf%n_ao)
!
      call output%print_matrix('debug', 'PCM Fock - Spin alpha + beta', &
                               wf%pcm_fock, wf%n_ao, wf%n_ao)
!
      call output%print_matrix('debug', 'QM/PCM Fock - Spin alpha', wf%ao_fock_a, wf%n_ao, wf%n_ao)
!
      call output%print_matrix('debug', 'QM/PCM Fock - Spin beta', wf%ao_fock_b, wf%n_ao, wf%n_ao)
!
!
   end subroutine update_fock_pcm_uhf
!
!
end submodule ao_fock
