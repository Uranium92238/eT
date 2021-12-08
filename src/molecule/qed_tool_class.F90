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
module qed_tool_class
!
!!
!!    QED tool class module
!!    Written by Tor S. Haugland, 2021
!!
!!    Class that stores information about the strong light-matter interaction.
!!    In particular, the coupling, polarization vectors, frequency and interaction
!!    integrals.
!
   use parameters
   use global_in,             only: input
   use global_out,            only: output
   use memory_manager_class,  only: mem
   use ao_tool_class,         only: ao_tool
   use array_utilities,       only: zero_array, symmetric_sandwich
!
   implicit none
!
   type :: qed_tool
!
      integer, private :: n_ao, n_el
      integer, public :: n_ph    ! n_photons per mode
      integer, public :: n_modes
!
      real(dp), dimension(:), allocatable, public :: frequency
      real(dp), dimension(:), allocatable, private :: coupling_bilinear, coupling_self ! n_modes
      real(dp), dimension(:, :), allocatable, private :: polarizations ! 3, n_modes
      real(dp), dimension(:), allocatable, private :: coherent_state ! n_modes
!
      real(dp), dimension(:,:,:), allocatable, private :: ao_eri_wx !needed for fast AO eri
      real(dp), dimension(:,:), allocatable, private :: ao_overlap
      real(dp), dimension(:,:,:), allocatable, private :: ao_mu
      real(dp), dimension(:,:,:), allocatable, private :: ao_q
      real(dp), dimension(:), allocatable, private :: mu_nuc
!
      logical, public :: is_optimizing_photons = .true. ! Optimize coherent state / photon displacement?
      logical, public :: is_basis_complete = .false. ! Approximate sum_p |p><p| -> 1? Removes references to virtual orbitals
!
   contains
!
      procedure, public :: initialize
      procedure, public :: cleanup
      procedure, public :: print_summary
!
      procedure, public :: add_ao_h_contribution
      procedure, public :: add_ao_g_contribution
      procedure, public :: get_ao_bilinear_integrals
      procedure, public :: update_coherent_state
!
      procedure, private :: get_polarized_ao_mu_wx
      procedure, private :: get_polarized_ao_q_wx
      procedure, private :: get_ao_self_energy_oei
      procedure, private :: initialize_ao_self_energy_eri
!
      procedure, private :: read_settings
      procedure, private, nopass :: wavevector_to_polarizations
!
   end type qed_tool
!
   interface qed_tool
!
      procedure :: new_qed_tool
!
   end interface qed_tool
!
!
contains
!
!
   function new_qed_tool() result(qed)
!!
!!    New QED tool
!!    Written by Tor S. Haugland, Oct 2021
!!
      implicit none
!
      type(qed_tool) :: qed
!
      call input%get_required_keyword('modes', 'qed', qed%n_modes)
!
      call mem%alloc(qed%frequency,             qed%n_modes)
      call mem%alloc(qed%polarizations,      3, qed%n_modes)
      call mem%alloc(qed%coherent_state,        qed%n_modes)
      call mem%alloc(qed%coupling_bilinear,     qed%n_modes)
      call mem%alloc(qed%coupling_self,         qed%n_modes)
!
      call qed%read_settings()
!
   end function new_qed_tool
!
!
   subroutine initialize(qed, ao)
!!
!!    Initialize
!!    Written by Tor S. Haugland, Dec 2021
!!
      use point_charges_class, only: point_charges
!
      implicit none
!
      class(qed_tool), intent(inout) :: qed
      class(ao_tool), intent(in) :: ao
!
      type(point_charges) :: pc
!
      qed%n_ao = ao%n
      qed%n_el = ao%get_n_electrons()
!
      call mem%alloc(qed%mu_nuc, 3)
      call ao%get_point_charges(pc)
      qed%mu_nuc = pc%get_dipole()
!
      call mem%alloc(qed%ao_overlap, qed%n_ao, qed%n_ao)
      call ao%get_oei('overlap', qed%ao_overlap)
!
      call mem%alloc(qed%ao_mu, qed%n_ao, qed%n_ao, 3)
      call ao%get_oei('dipole', qed%ao_mu)
!
      if (qed%is_basis_complete) then
         call mem%alloc(qed%ao_q, qed%n_ao, qed%n_ao, 6)
         call ao%get_oei('quadrupole', qed%ao_q)
      endif
!
      call mem%alloc(qed%ao_eri_wx, qed%n_ao, qed%n_ao,  qed%n_modes)
      call qed%initialize_ao_self_energy_eri()
!
   end subroutine initialize
!
!
   subroutine cleanup(qed)
!!
!!    Cleanup
!!    Written by Tor S. Haugland, Oct 2021
!!
      implicit none
!
      class(qed_tool), intent(inout) :: qed
!
      call mem%dealloc(qed%frequency,         qed%n_modes)
      call mem%dealloc(qed%polarizations,  3, qed%n_modes)
      call mem%dealloc(qed%coherent_state,    qed%n_modes)
      call mem%dealloc(qed%coupling_bilinear, qed%n_modes)
      call mem%dealloc(qed%coupling_self,     qed%n_modes)
      call mem%dealloc(qed%ao_eri_wx, qed%n_ao, qed%n_ao, qed%n_modes)
!
      call mem%dealloc(qed%mu_nuc, 3)
      call mem%dealloc(qed%ao_overlap, qed%n_ao, qed%n_ao)
      call mem%dealloc(qed%ao_mu, qed%n_ao, qed%n_ao, 3)
      if (qed%is_basis_complete) &
            call mem%dealloc(qed%ao_q, qed%n_ao, qed%n_ao, 6)
!
   end subroutine cleanup
!
!
   subroutine read_settings(qed)
!!
!!    Read settings
!!    Written by Tor S. Haugland, Sep 2018
!!
!!    Designed to be overwritten by descendants.
!!
      implicit none
!
      class(qed_tool), intent(inout) :: qed
      real(dp), dimension(3) :: temp
      real(dp), dimension(:), allocatable :: coupling
      integer :: k
!
      if (input%is_keyword_present('frequency', 'qed')) then
         call input%get_array_for_keyword('frequency', 'qed', qed%n_modes, qed%frequency)
      else
         call output%error_msg("Keyword frequency must be specified.")
      endif
!
      if (input%is_keyword_present('wavevector', 'qed')) then
         call input%get_array_for_keyword('wavevector', 'qed', 3, temp)
         do k = 1, qed%n_modes, 2
            call qed%wavevector_to_polarizations(temp, qed%polarizations(:,k), qed%polarizations(:,k+1))
         enddo
      else if (input%is_keyword_present('polarization', 'qed')) then
         call input%get_array_for_keyword('polarization', 'qed', 3, temp)
         do k = 1, qed%n_modes
            qed%polarizations(:,k) = temp(:)
         enddo
      else
         call output%error_msg("Keyword wavevector or polarization must be specified")
      endif
!
      qed%coherent_state(:) = 0.0d0
      if ( input%is_keyword_present('coherent state', 'qed') &
           .or. input%is_keyword_present('hf coherent state', 'qed')) then
         qed%is_optimizing_photons = .false.
         call input%get_array_for_keyword('coherent state', 'qed', qed%n_modes, qed%coherent_state)
         call input%get_array_for_keyword('hf coherent state', 'qed', qed%n_modes, qed%coherent_state)
      endif
!
!     Either coupling or bilinear and self couplings must be specified
!
      call mem%alloc(coupling, qed%n_modes)
      call zero_array(coupling, qed%n_modes)
!
      if ( input%is_keyword_present('coupling', 'qed')) then
         call input%get_array_for_keyword('coupling', 'qed', qed%n_modes, coupling)
      else
         call output%warning_msg("Coupling is not specified and is set to zero.")
      endif
!
      do k = 1, qed%n_modes
         qed%coupling_bilinear(k) = coupling(k) * sqrt(qed%frequency(k) * half)
         qed%coupling_self(k) = coupling(k)**2 * half
      enddo
      call mem%dealloc(coupling, qed%n_modes)
!
      if ( input%is_keyword_present('coupling bilinear', 'qed')) &
         call input%get_array_for_keyword('coupling bilinear', 'qed', qed%n_modes, qed%coupling_bilinear)
      if ( input%is_keyword_present('coupling self', 'qed')) &
         call input%get_array_for_keyword('coupling self', 'qed', qed%n_modes, qed%coupling_self)
!
      qed%is_basis_complete = input%is_keyword_present('quadrupole oei', 'qed')
!
   end subroutine read_settings
!
!
   subroutine initialize_ao_self_energy_eri(qed)
!!
!!    Initialize AO self-energy ERI
!!    Written by T. S. Haugland, Oct 2021
!!
!!    Sets the dipole integrals needed for fast computation
!!    of two electron integrals (ERI). Computing g_wx for every get_g
!!    call is too expensive.
!!
      implicit none
!
      class(qed_tool), intent(inout) :: qed
      integer :: mode
!
      do mode = 1, qed%n_modes
         call qed%get_polarized_ao_mu_wx(qed%ao_eri_wx(:,:,mode), mode)
         call dscal(qed%n_ao**2, sqrt(two*qed%coupling_self(mode)), qed%ao_eri_wx(:,:,mode), 1)
      enddo
!
   end subroutine initialize_ao_self_energy_eri
!
!
   subroutine get_polarized_ao_mu_wx(qed, mu_wx, mode)
!!
!!    Get polarized AO mu_wx
!!    Written by T. S. Haugland, Sep 2018
!!
!!    Gets the product between the photon transversal polarization
!!    and the dipole moments.
!!       mu_wx = <w| (e dot mu_el) |x> + (e dot mu_nuc) S_wx / N_el
!!          e : transversal polarization vector
!!
      implicit none
!
      class(qed_tool), intent(in) :: qed
      real(dp), dimension(:,:), intent(out) :: mu_wx
      integer, intent(in) :: mode
!
      real(dp) :: d_nuc
      real(dp), external :: ddot
      integer :: i
!
!     Electronic Contribution
!
      call zero_array(mu_wx, qed%n_ao**2)
      do i = 1, 3
         call daxpy(qed%n_ao**2, qed%polarizations(i,mode), qed%ao_mu(:,:,i), 1, mu_wx, 1)
      enddo
!
!     Nuclear Contribution
!
      d_nuc = ddot(3, qed%mu_nuc, 1, qed%polarizations(:,mode), 1)
      call daxpy(qed%n_ao**2, d_nuc / qed%n_el, qed%ao_overlap, 1, mu_wx, 1)
!
   end subroutine get_polarized_ao_mu_wx
!
!
   subroutine get_polarized_ao_q_wx(qed, Q_wx, mode)
!!
!!    Get polarized AO Q_wx
!!    Written by T. S. Haugland, Nov 2021
!!
!!    Gets the product between the photon transversal polarization
!!    and the quadrupole moments. Includes nuclear contribution.
!!
!!    Q_wx = <w| (e^T dot Q dot e) |x>
!!       e : transversal polarization vector
!!
      implicit none
!
      class(qed_tool), intent(in) :: qed
      real(dp), dimension(:,:), intent(out) :: Q_wx
      integer, intent(in) :: mode
!
      real(dp) :: d_nuc
      real(dp), external :: ddot
      integer :: i
!
      call zero_array(Q_wx, qed%n_ao**2)
      call daxpy(qed%n_ao**2, -qed%polarizations(1,mode) * qed%polarizations(1,mode),       qed%ao_q(:,:,1), 1, Q_wx, 1) !xx
      call daxpy(qed%n_ao**2, -qed%polarizations(1,mode) * qed%polarizations(2,mode) * two, qed%ao_q(:,:,2), 1, Q_wx, 1) !xy
      call daxpy(qed%n_ao**2, -qed%polarizations(1,mode) * qed%polarizations(3,mode) * two, qed%ao_q(:,:,3), 1, Q_wx, 1) !xz
      call daxpy(qed%n_ao**2, -qed%polarizations(2,mode) * qed%polarizations(2,mode),       qed%ao_q(:,:,4), 1, Q_wx, 1) !yy
      call daxpy(qed%n_ao**2, -qed%polarizations(2,mode) * qed%polarizations(3,mode) * two, qed%ao_q(:,:,5), 1, Q_wx, 1) !yz
      call daxpy(qed%n_ao**2, -qed%polarizations(3,mode) * qed%polarizations(3,mode),       qed%ao_q(:,:,6), 1, Q_wx, 1) !zz
!
      d_nuc = ddot(3, qed%mu_nuc, 1, qed%polarizations(:,mode), 1)
      call daxpy(qed%n_ao**2, (d_nuc/qed%n_el)**2, qed%ao_overlap, 1, Q_wx, 1)
!
      do i = 1, 3
         call daxpy(qed%n_ao**2, 2*(d_nuc/qed%n_el)*qed%polarizations(i,mode), qed%ao_mu(:,:,i), 1, Q_wx, 1)
      enddo
!
   end subroutine get_polarized_ao_q_wx
!
!
   subroutine get_ao_bilinear_integrals(qed, g_wx, mode)
!!
!!    Get AO bilinear integrals
!!    Written by T. S. Haugland, June 2021
!!
!!    Gets the bilinear interaction term in the AO basis,
!!    the g in g(b+b^+).
!!       g_wx = mu_wx * lambda * sqrt(omega/2)
!!
      implicit none
!
      class(qed_tool), intent(in) :: qed
      real(dp), dimension(:,:), intent(out) :: g_wx
      integer, intent(in) :: mode
!
      call qed%get_polarized_ao_mu_wx(g_wx, mode)
      call dscal(qed%n_ao**2, qed%coupling_bilinear(mode), g_wx, 1)
!
   end subroutine get_ao_bilinear_integrals
!
!
   subroutine get_ao_self_energy_oei(qed, g_wx)
!!
!!    Get AO self-energy OEI
!!    Written by T. S. Haugland, June 2021
!!
!!    Gets the quadratic one-electron interaction in the AO basis,
!!    also known as the dipole self-energy.
!!       g2_wx = sum_p g_wp g_px = lambda^2 / 2 * sum_p mu_wp mu_px
!!
      use array_utilities, only: invert
!
      implicit none
!
      class(qed_tool), intent(in) :: qed
      real(dp), dimension(:,:), intent(out) :: g_wx
!
      real(dp), dimension(:,:), allocatable :: D_full, mu_wx, Q_wx, muDmu_wx
      integer :: mode
!
      if (qed%is_basis_complete) then
!
!        Assume sum_p |p><p| = 1 to compute quadrupole integrals instead of full density matrix
!
         call mem%alloc(Q_wx, qed%n_ao, qed%n_ao)
!
         call zero_array(g_wx, qed%n_ao**2)
         do mode = 1, qed%n_modes
            call qed%get_polarized_ao_q_wx(Q_wx, mode)
            call daxpy(qed%n_ao**2, qed%coupling_self(mode), Q_wx, 1, g_wx, 1)
         enddo
!
         call mem%dealloc(Q_wx, qed%n_ao, qed%n_ao)
!
      else
!
!        Invert S to compute the full density matrix
!
         call mem%alloc(D_full, qed%n_ao, qed%n_ao)
         call invert(D_full, qed%ao_overlap, qed%n_ao)
!
         call mem%alloc(mu_wx, qed%n_ao, qed%n_ao)
         call mem%alloc(muDmu_wx, qed%n_ao, qed%n_ao)
!
         call zero_array(g_wx, qed%n_ao**2)
         do mode = 1, qed%n_modes
            call qed%get_polarized_ao_mu_wx(mu_wx, mode)
            call symmetric_sandwich(muDmu_wx, D_full, mu_wx, qed%n_ao, qed%n_ao)
            call daxpy(qed%n_ao**2, qed%coupling_self(mode), muDmu_wx, 1, g_wx, 1)
         enddo
!
         call mem%dealloc(D_full, qed%n_ao, qed%n_ao)
         call mem%dealloc(mu_wx, qed%n_ao, qed%n_ao)
         call mem%dealloc(muDmu_wx, qed%n_ao, qed%n_ao)
!
      endif
!
   end subroutine get_ao_self_energy_oei
!
!
   subroutine add_ao_h_contribution(qed, h_wx)
!!
!!    Add AO h contribution
!!    Written by T. S. Haugland, Oct 2021
!!
!!    Adds the QED AO contribution to h.
!!       h_wx += [g^2]_wx / omega + 2z * g_wx + omega z^2 S_wx / N_e
!!
      implicit none
!
      class(qed_tool), intent(in) :: qed
      real(dp), dimension(:,:), intent(out) :: h_wx
      real(dp), dimension(:,:), allocatable :: g_wx
      integer :: mode
!
      call mem%alloc(g_wx, qed%n_ao, qed%n_ao)
      call qed%get_ao_self_energy_oei(g_wx)
      call daxpy(qed%n_ao**2, one, g_wx, 1, h_wx, 1)
!
      do mode = 1, qed%n_modes
         call qed%get_ao_bilinear_integrals(g_wx, mode)
         call daxpy(qed%n_ao**2, two * qed%coherent_state(mode), g_wx, 1, h_wx, 1)
      enddo
      call mem%dealloc(g_wx, qed%n_ao, qed%n_ao)
!
      call daxpy(qed%n_ao**2, sum(qed%frequency * qed%coherent_state**2)/ qed%n_el, &
         qed%ao_overlap, 1, h_wx, 1)
!
   end subroutine add_ao_h_contribution
!
!
   subroutine add_ao_g_contribution(qed, g, shellA, shellB, shellC, shellD)
!!
!!    Add AO g contribution
!!    Written by Tor S. Haugland, Oct 2021
!!
!!    Add quadratic self-energy QED contribution to electron repulsion integrals for the
!!    shell quartet (A, B, C, D).
!!
!!       g += 2/omega * g_wx g_yz
!!
      use named_range_class
!
      implicit none
!
      class(qed_tool), intent(in) :: qed
      type(range_), intent(in) :: shellA, shellB, shellC, shellD
      real(dp), dimension(shellA%length * shellB%length * &
                          shellC%length * shellD%length), target, intent(out) :: g
!
      real(dp), dimension(:,:,:,:), pointer :: g_p
      integer :: w, x, y, z, mode
!
      g_p(shellA%first : shellA%get_last(), &
          shellB%first : shellB%get_last(), &
          shellC%first : shellC%get_last(), &
          shellD%first : shellD%get_last()) => g
!
      do mode = 1, qed%n_modes
         do z = shellD%first, shellD%get_last()
            do y = shellC%first, shellC%get_last()
               do x = shellB%first, shellB%get_last()
                  do w = shellA%first, shellA%get_last()
                     g_p(w,x,y,z) = g_p(w,x,y,z) + qed%ao_eri_wx(w,x,mode) * qed%ao_eri_wx(y,z,mode)
                  enddo
               enddo
            enddo
         enddo
      enddo
!
   end subroutine add_ao_g_contribution
!
!
   subroutine print_summary(qed)
!!
!!    Print summary
!!    Written by T. S. Haugland, Oct 2021
!!
      implicit none
!
      class(qed_tool), intent(in) :: qed
      integer :: mode
!
      call output%printf('m', '- QED Parameters and properties', ffs='(/t3,a)')
!
      call output%printf('m', 'Optimize photons?   :    (l0)', logs=[qed%is_optimizing_photons], ffs='(/t6,a)')
      call output%printf('m', 'Complete basis?     :    (l0)', logs=[qed%is_basis_complete], ffs='(t6,a)')
!
      do mode = 1, qed%n_modes
         call output%printf('m', 'Mode (i0)', ints=[mode], ffs='(/t6,a)')
         call output%printf('m', ' Frequency          : (f17.12)', reals=[qed%frequency(mode)], ffs='(t6,a)')
         call output%printf('m', ' Polarization       : (f17.12) (f17.12) (f17.12)', reals=qed%polarizations(:,mode), ffs='(t6,a)')
         call output%printf('m', ' Coupling', ffs='(t6,a)')
         call output%printf('m', '  Bilinear          : (f17.12)', reals=[qed%coupling_bilinear(mode)], ffs='(t6,a)')
         call output%printf('m', '  Quadratic         : (f17.12)', reals=[qed%coupling_self(mode)], ffs='(t6,a)')
         call output%printf('m', ' Coherent state     : (f17.12)', reals=[qed%coherent_state(mode)], ffs='(t6,a)')
      enddo
!
   end subroutine print_summary
!
!
   subroutine wavevector_to_polarizations(wavevector, pol_1, pol_2)
!!
!!    Wavevector to polarizations
!!    Written by Tor S. Haugland, May 2021
!!
!!    From a wavevector, form two orthogonal polarization directions.
!!
      use math_utilities, only: cross_product_R3
!
      implicit none
!
      real(dp), dimension(3), intent(in) :: wavevector
      real(dp), dimension(3), intent(out) :: pol_1, pol_2
      real(dp), dimension(3) :: rand_vec
      logical :: is_parallel
!
      rand_vec = (/ 0.0d0, 0.0d0, 1.0d0 /)
      is_parallel = norm2(abs(rand_vec)-abs(wavevector)) < 1d-10
      if (is_parallel) rand_vec = (/ 0.0d0, 1.0d0, 0.0d0 /)
!
      pol_1 = cross_product_R3(rand_vec, wavevector)
      pol_1 = pol_1 / norm2(pol_1)
!
      pol_2 = cross_product_R3(wavevector, pol_1)
      pol_2 = pol_2 / norm2(pol_2)
!
   end subroutine wavevector_to_polarizations
!
!
   subroutine update_coherent_state(qed, ao_density)
!!
!!    Update coherent state
!!    Written by T. S. Haugland, June 2021
!!
!!    Coherent state z is given from the bilinear coupling g_b and frequency w,
!!       z = - <g_b> / w = - Tr[D g_b] / w
!!
      implicit none
!
      class(qed_tool) :: qed
      real(dp), dimension(qed%n_ao, qed%n_ao) :: ao_density

      real(dp), external :: ddot
      real(dp), dimension(:,:), allocatable :: g_wx
      integer :: mode
!
      call mem%alloc(g_wx, qed%n_ao, qed%n_ao)
      do mode = 1, qed%n_modes
         call qed%get_ao_bilinear_integrals(g_wx, mode)
         qed%coherent_state(mode) = - ddot(qed%n_ao**2, ao_density, 1, g_wx, 1) / qed%frequency(mode)
      enddo
      call mem%dealloc(g_wx, qed%n_ao, qed%n_ao)
!
   end subroutine update_coherent_state
!
!
end module qed_tool_class
