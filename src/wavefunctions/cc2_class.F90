module cc2_class
!
!!
!!    Coupled cluster singles and perturbative doubles (CC2) class module
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto and Alexander Paul, 2018
!!
!
   use ccs_class
!
   implicit none
!
   type, extends(ccs) :: cc2
!
      real(dp), dimension(:,:,:,:), allocatable :: u
!
   contains
!
      procedure :: prepare                         => prepare_cc2
!
      procedure :: construct_u                     => construct_u_cc2
!
      procedure :: construct_omega                 => construct_omega_cc2
!
      procedure :: omega_cc2_a1                    => omega_cc2_a1_cc2
      procedure :: omega_cc2_b1                    => omega_cc2_b1_cc2
      procedure :: omega_cc2_c1                    => omega_cc2_c1_cc2
!
      procedure :: calculate_energy                => calculate_energy_cc2
!
      procedure :: prepare_for_jacobian            => prepare_for_jacobian_cc2
!
      procedure :: jacobian_transform_trial_vector => jacobian_transform_trial_vector_cc2
      procedure :: jacobian_cc2_transformation     => jacobian_cc2_transformation_cc2
!
      procedure :: jacobian_cc2_a1                 => jacobian_cc2_a1_cc2
      procedure :: jacobian_cc2_b1                 => jacobian_cc2_b1_cc2
      procedure :: jacobian_cc2_a2                 => jacobian_cc2_a2_cc2
      procedure :: jacobian_cc2_b2                 => jacobian_cc2_b2_cc2
!
      procedure :: initialize_u                    => initialize_u_cc2 
      procedure :: destruct_u                      => destruct_u_cc2 
!
      procedure :: initialize_amplitudes           => initialize_amplitudes_cc2 
      procedure :: destruct_amplitudes             => destruct_amplitudes_cc2 
!
   end type cc2
!
   interface
!
      include "../submodules/cc2/omega_cc2_interface.F90"
      include "../submodules/cc2/jacobian_cc2_interface.F90"
!
   end interface 
!
!
contains
!
!
   subroutine prepare_cc2(wf, ref_wf)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(cc2) :: wf
!
      class(hf) :: ref_wf
!
      integer(i15) :: p
!
      wf%name = 'cc2'
!
      wf%system = ref_wf%system
!
      wf%n_ao   = ref_wf%n_ao
      wf%n_mo   = ref_wf%n_mo
      wf%n_o    = ref_wf%n_o
      wf%n_v    = ref_wf%n_v
!
      wf%hf_energy = ref_wf%energy
!
      wf%n_t1         = (wf%n_o)*(wf%n_v)
      wf%n_amplitudes = wf%n_t1
!
      call wf%initialize_fock_ij()
      call wf%initialize_fock_ia()
      call wf%initialize_fock_ai()
      call wf%initialize_fock_ab()
!
      call wf%initialize_fock_diagonal()
!
      wf%fock_ij(:,:) = ref_wf%mo_fock(1 : wf%n_o, 1 : wf%n_o)
      wf%fock_ia(:,:) = ref_wf%mo_fock(1 : wf%n_o, wf%n_o + 1 : wf%n_mo)
      wf%fock_ai(:,:) = ref_wf%mo_fock(wf%n_o + 1 : wf%n_mo, 1 : wf%n_o)
      wf%fock_ab(:,:) = ref_wf%mo_fock(wf%n_o + 1 : wf%n_mo, wf%n_o + 1 : wf%n_mo)
!
      do p = 1, wf%n_mo
!
         wf%fock_diagonal(p, 1) = ref_wf%mo_fock(p, p)
!
      enddo
!
      call ref_wf%mo_transform_and_save_h()
!
      call wf%initialize_orbital_coefficients()
      wf%orbital_coefficients = ref_wf%orbital_coefficients
!
   end subroutine prepare_cc2
!
!
   module subroutine calculate_energy_cc2(wf)
!!
!!    Calculate energy 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
!!    E = E_HF + sum_aibj (t_i^a*t_j^b + t_ij^ab) L_iajb
!!
      class(cc2), intent(inout) :: wf 
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj, g_iajb 
!
      real(dp) :: correlation_energy
!
      integer(i15) :: a, i, b, j
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_vovo(g_aibj)
      call wf%get_ovov(g_iajb)
!
      correlation_energy = zero 
!
!$omp parallel do private(a,i,b,j) reduction(+:correlation_energy)
      do b = 1, wf%n_v
         do i = 1, wf%n_o 
            do j = 1, wf%n_o 
               do a = 1, wf%n_v
!
                  correlation_energy = correlation_energy +                                &
                                       (wf%t1(a, i)*wf%t1(b, j) -                          &
                                       (g_aibj(a,i,b,j))/(wf%fock_diagonal(wf%n_o + a, 1)  &
                                                         + wf%fock_diagonal(wf%n_o + b, 1) &
                                                         - wf%fock_diagonal(i,1)           &
                                                         - wf%fock_diagonal(j,1)))         &
                                       *(two*g_iajb(i,a,j,b)-g_iajb(i,b,j,a))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do 
!
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      wf%energy = wf%hf_energy + correlation_energy
!
   end subroutine calculate_energy_cc2
!
!
   subroutine construct_u_cc2(wf)
!!
!!    Construct U 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
!!    Construct
!!
!!       u_aibj = 2t_aibj - t_ajbi
!!
!!    with
!!
!!       t_aibj = - g_aibj/ε_aibj
!!
!!    where
!!
!!       ε_aibj = ε_a - ε_i + ε_b - ε_j 
!!
!!    and ε_r is the r'th orbital energy.
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
!
      integer(i15) :: a, i, b, j
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)  
!
      call wf%get_vovo(g_aibj)
!
      do b = 1, wf%n_v 
         do j = 1, wf%n_o 
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  wf%u(a, i, b, j) = (two*g_aibj(a, i, b, j) - g_aibj(a, j, b, i))/ &
                                          (wf%fock_diagonal(i, 1) + wf%fock_diagonal(j, 1) &
                                        - wf%fock_diagonal(wf%n_o + a, 1) &
                                        - wf%fock_diagonal(wf%n_o + b, 1) )
!
               enddo
            enddo
         enddo 
      enddo
!    
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)      
!
   end subroutine construct_u_cc2
!
!
   subroutine initialize_u_cc2(wf)
!!
!!    Initialize u 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      if (.not. allocated(wf%u)) call mem%alloc(wf%u, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine initialize_u_cc2
!
!
   subroutine destruct_u_cc2(wf)
!!
!!    Initialize u 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      if (allocated(wf%u)) call mem%dealloc(wf%u, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine destruct_u_cc2
!
!
   subroutine initialize_amplitudes_cc2(wf)
!!
!!    Initialize amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      call wf%initialize_t1()
      call wf%initialize_u()
!
   end subroutine initialize_amplitudes_cc2
!
!
   subroutine destruct_amplitudes_cc2(wf)
!!
!!    Destruct amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      call wf%destruct_t1()
      call wf%destruct_u()
!
   end subroutine destruct_amplitudes_cc2
!
end module cc2_class
