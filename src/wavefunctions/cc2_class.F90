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
      integer :: n_t2
!
      real(dp), dimension(:,:,:,:), allocatable :: u
!
   contains
!
      procedure :: prepare                                     => prepare_cc2
!
      procedure :: construct_u                                 => construct_u_cc2
!
      procedure :: construct_omega                             => construct_omega_cc2
!
      procedure :: omega_cc2_a1                                => omega_cc2_a1_cc2
      procedure :: omega_cc2_b1                                => omega_cc2_b1_cc2
      procedure :: omega_cc2_c1                                => omega_cc2_c1_cc2
!
      procedure :: calculate_energy                            => calculate_energy_cc2
!
      procedure :: prepare_for_jacobian                        => prepare_for_jacobian_cc2
!
      procedure :: jacobian_transform_trial_vector             => jacobian_transform_trial_vector_cc2
      procedure :: jacobian_cc2_transformation                 => jacobian_cc2_transformation_cc2
!
      procedure :: jacobian_cc2_a1                             => jacobian_cc2_a1_cc2
      procedure :: jacobian_cc2_b1                             => jacobian_cc2_b1_cc2
      procedure :: jacobian_cc2_a2                             => jacobian_cc2_a2_cc2
      procedure :: jacobian_cc2_b2                             => jacobian_cc2_b2_cc2
!
      procedure :: prepare_for_jacobian_transpose              => prepare_for_jacobian_transpose_cc2
!
      procedure :: jacobian_transpose_transform_trial_vector   => jacobian_transpose_transform_trial_vector_cc2
      procedure :: jacobian_transpose_cc2_transformation       => jacobian_transpose_cc2_transformation_cc2
!
      procedure :: jacobian_transpose_cc2_a1                   => jacobian_transpose_cc2_a1_cc2
      procedure :: jacobian_transpose_cc2_b1                   => jacobian_transpose_cc2_b1_cc2
      procedure :: jacobian_transpose_cc2_a2                   => jacobian_transpose_cc2_a2_cc2
      procedure :: jacobian_transpose_cc2_b2                   => jacobian_transpose_cc2_b2_cc2
!
      procedure :: initialize_u                                => initialize_u_cc2 
      procedure :: destruct_u                                  => destruct_u_cc2 
!
      procedure :: initialize_amplitudes                       => initialize_amplitudes_cc2 
      procedure :: destruct_amplitudes                         => destruct_amplitudes_cc2 
!
      procedure :: get_es_orbital_differences                  => get_es_orbital_differences_cc2
!
      procedure :: construct_multiplier_equation               => construct_multiplier_equation_cc2
!
      procedure :: get_cvs_projector                           => get_cvs_projector_cc2
!
   end type cc2
!
   interface
!
      include "../submodules/cc2/omega_cc2_interface.F90"
      include "../submodules/cc2/jacobian_cc2_interface.F90"
      include "../submodules/cc2/jacobian_transpose_cc2_interface.F90"
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
      integer :: p
!
      wf%name_ = 'cc2'
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
      wf%n_t1            = (wf%n_o)*(wf%n_v)
      wf%n_t2            = wf%n_t1*(wf%n_t1+1)/2
      wf%n_gs_amplitudes = wf%n_t1
      wf%n_es_amplitudes = wf%n_t1 + wf%n_t2
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
      integer :: a, i, b, j
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
      integer :: a, i, b, j
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
!
   subroutine get_es_orbital_differences_cc2(wf, orbital_differences, N)
!!
!!    Get orbital differences 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      integer, intent(in) :: N 
      real(dp), dimension(N), intent(inout) :: orbital_differences
!
      integer :: a, i, ai, b, j, bj, aibj
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) 
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            orbital_differences(ai) = wf%fock_diagonal(a + wf%n_o, 1) - wf%fock_diagonal(i, 1)
!
            do j = 1, wf%n_o 
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j-1) + b 
!
                  if (ai .ge. bj) then
!
                     aibj = (ai*(ai-3)/2) + ai + bj
!
                     orbital_differences(aibj + (wf%n_o)*(wf%n_v)) = wf%fock_diagonal(a + wf%n_o, 1) - wf%fock_diagonal(i, 1) &
                                                                      +  wf%fock_diagonal(b + wf%n_o, 1) - wf%fock_diagonal(j, 1)
!
                  endif
!
               enddo
            enddo  
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine get_es_orbital_differences_cc2
!
!
   subroutine construct_multiplier_equation_cc2(wf, equation)
!!
!!    Construct multiplier equation 
!!    Written by Sarai D. Folkestad, Feb 2019
!!
!!    Constructs 
!!
!!       t-bar^T A + eta,
!!
!!    and places the result in 'equation'.
!!
!!    Solves analytically for tbar_aibj
!!
!!       tbar_aibj = - (η_aibj + sum_ai tbar_ai A_ai,aibj)/ε_aibj
!!
!!    where
!!
!!       η_aibj = 2 L_iajb       
!!
!!    and uses this to set up 'equation'
!!
!!       η_ai + sum_bj tbar_bj A_bj,ai + sum_bjck tbar_bjck A_{bjck,ai}
!!
      implicit none 
!
      class(cc2), intent(in) :: wf 
!
      real(dp), dimension(wf%n_gs_amplitudes, 1), intent(inout) :: equation 
!
      real(dp), dimension(:,:), allocatable :: eta 
      real(dp), dimension(:,:,:,:), allocatable :: t2bar
      real(dp), dimension(:,:,:,:), allocatable :: g_iajb
!
      integer :: a, b, i, j
!
!     Construct t2bar
!
      call mem%alloc(t2bar, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      t2bar = zero
!
!     t2bar = sum_ai tbar_ai A_ai,aibj
!
      call wf%jacobian_transpose_cc2_a2(t2bar, wf%t1bar)
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%get_ovov(g_iajb)
!
!     t2bar += η_aibj
!
      call add_2143_to_1234(four, g_iajb, t2bar, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2341_to_1234(-two, g_iajb, t2bar, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     t2bar = t2bar/(-ε_aibj)
!
!$omp parallel do private(a, b, i, j)
      do b = 1, wf%n_v
         do j = 1, wf%n_o
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  t2bar(a, i, b, j) = t2bar(a, i, b, j)/(- wf%fock_diagonal(a + wf%n_o, 1) &
                                                        -  wf%fock_diagonal(b + wf%n_o, 1) &
                                                        +  wf%fock_diagonal(i, 1) &
                                                        +  wf%fock_diagonal(j, 1))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Set up the multipliers equation
!
      equation = zero
!
!     equation += sum_bjck tbar_bjck A_{bjck,ai}
!
      call wf%jacobian_transpose_cc2_b1(equation, t2bar)
!
      call mem%dealloc(t2bar, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     equation += sum_bj tbar_bj A_bj,ai
!  
      call wf%jacobian_transpose_ccs_a1(equation, wf%t1bar)
      call wf%jacobian_transpose_ccs_b1(equation, wf%t1bar)
      call wf%jacobian_transpose_cc2_a1(equation, wf%t1bar)
!
!     Add eta, equation = t-bar^T A + eta 
!
      call mem%alloc(eta, wf%n_gs_amplitudes, 1)
      call wf%construct_eta(eta)
!
      call daxpy(wf%n_gs_amplitudes, one, eta, 1, equation, 1)
!
      call mem%dealloc(eta, wf%n_gs_amplitudes, 1)
!
   end subroutine construct_multiplier_equation_cc2
!
!
   subroutine get_cvs_projector_cc2(wf, projector, n_cores, core_MOs)
!!
!!    Get CVS projector
!!    Written by Sarai D. Folekstad, Oct 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(out) :: projector
!
      integer, intent(in) :: n_cores
!
      integer, dimension(n_cores, 1), intent(in) :: core_MOs
!
      integer :: core, i, a, ai, j, b, bj, aibj
!
      projector = zero
!
      do core = 1, n_cores
!
        i = core_MOs(core, 1)
!
!$omp parallel do private (a, ai, j, b, bj, aibj)
        do a = 1, wf%n_v
!
           ai = wf%n_v*(i - 1) + a
           projector(ai, 1) = one
!
            do j = 1, wf%n_o 
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
                  aibj = max(ai, bj)*(max(ai, bj) - 3)/2 + ai + bj
!                  
                  projector(aibj + (wf%n_o)*(wf%n_v), 1) = one
!
               enddo
            enddo
        enddo
!$omp end parallel do
!
     enddo
!
   end subroutine get_cvs_projector_cc2
!
!
end module cc2_class
