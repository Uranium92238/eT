module lowmem_cc2_class
!
!!
!!    Low-memory coupled cluster singles and perturbative doubles (CC2) class module
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto and Alexander Paul, 2018
!!
!
   use ccs_class
!
   implicit none
!
   type, extends(ccs) :: lowmem_cc2
!
   contains
!
      procedure :: prepare          => prepare_lowmem_cc2
!
      procedure :: construct_omega  => construct_omega_lowmem_cc2
!
      procedure :: omega_cc2_a1     => omega_cc2_a1_lowmem_cc2
      procedure :: omega_cc2_b1     => omega_cc2_b1_lowmem_cc2
      procedure :: omega_cc2_c1     => omega_cc2_c1_lowmem_cc2
!
      procedure :: calculate_energy => calculate_energy_lowmem_cc2
!
      procedure :: construct_excited_state_equation => construct_excited_state_equation_lowmem_cc2
!
      procedure :: effective_jacobian_transformation => effective_jacobian_transformation_lowmem_cc2
!
      procedure :: jacobian_cc2_a1 => jacobian_cc2_a1_lowmem_cc2
      procedure :: jacobian_cc2_b1 => jacobian_cc2_b1_lowmem_cc2
!
      procedure :: effective_jacobian_cc2_a1 => effective_jacobian_cc2_a1_lowmem_cc2
      procedure :: effective_jacobian_cc2_b1 => effective_jacobian_cc2_b1_lowmem_cc2
      procedure :: effective_jacobian_cc2_c1 => effective_jacobian_cc2_c1_lowmem_cc2
      procedure :: effective_jacobian_cc2_d1 => effective_jacobian_cc2_d1_lowmem_cc2
      procedure :: effective_jacobian_cc2_e1 => effective_jacobian_cc2_e1_lowmem_cc2
      procedure :: effective_jacobian_cc2_f1 => effective_jacobian_cc2_f1_lowmem_cc2
!
   end type lowmem_cc2
!
   interface
!
      include "../submodules/lowmem_cc2/omega_lowmem_cc2_interface.F90"
      include "../submodules/lowmem_cc2/jacobian_lowmem_cc2_interface.F90"
!
   end interface
!
!
contains
!
!
   subroutine prepare_lowmem_cc2(wf, ref_wf)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(lowmem_cc2) :: wf
!
      class(hf) :: ref_wf
!
      integer :: p
!
      wf%name_ = 'low memory cc2'
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
      wf%n_gs_amplitudes = wf%n_t1
      wf%n_es_amplitudes = wf%n_t1
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
   end subroutine prepare_lowmem_cc2
!
!
   subroutine calculate_energy_lowmem_cc2(wf)
!!
!!     Calculate energy (lowmem_CC2)
!!     Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!     Andreas Skeidsvoll, 2018
!!
!!     Calculates the lowmem_CC2 energy. This is only equal to the actual
!!     energy when the ground state equations are solved, of course.
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
      real(dp), dimension(:,:,:,:), allocatable :: g_iajb
!
      real(dp) :: correlation_energy
!
      integer :: i, j, a, b
!
      integer :: req0, req1_i, req1_j, req2
!
      integer :: current_i_batch, current_j_batch
!
      type(batching_index) :: batch_i, batch_j
!
      req0 = 0
!
      req1_i = (wf%n_v)*(wf%integrals%n_J)
      req1_j = (wf%n_v)*(wf%integrals%n_J)
!
      req2 =  2*(wf%n_v**2)
!
      call batch_i%init(wf%n_o)
      call batch_j%init(wf%n_o)
!
      call mem%batch_setup(batch_i, batch_j, req0, req1_i, req1_j, req2)
!
      correlation_energy = zero
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(current_j_batch)
!
            call mem%alloc(g_aibj, wf%n_v, batch_i%length, wf%n_v, batch_j%length)
            call mem%alloc(g_iajb, batch_i%length, wf%n_v, batch_j%length, wf%n_v)
!
            call wf%get_vovo(g_aibj, &
                              1, wf%n_v, &
                              batch_i%first, batch_i%last, &
                              1, wf%n_v, &
                              batch_j%first, batch_j%last)
!
            call wf%get_ovov(g_iajb, &
                              batch_i%first, batch_i%last, &
                              1, wf%n_v, &
                              batch_j%first, batch_j%last, &
                              1, wf%n_v)
!
!$omp parallel do private(b,i,j,a) reduction(+:correlation_energy)
            do b = 1, wf%n_v
               do i = 1, batch_i%length
                  do j = 1, batch_j%length
                     do a = 1, wf%n_v
!
                        correlation_energy = correlation_energy + &
                                             (wf%t1(a, i + batch_i%first - 1)*wf%t1(b, j + batch_j%first - 1) &
                                             - (g_aibj(a, i, b, j))/(wf%fock_diagonal(wf%n_o + a, 1) &
                                                               + wf%fock_diagonal(wf%n_o + b, 1) &
                                                               - wf%fock_diagonal(i + batch_i%first - 1,1) &
                                                               - wf%fock_diagonal(j + batch_j%first - 1,1)))&
                                             *(two*g_iajb(i, a, j, b)-g_iajb(i, b, j, a))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(g_aibj, wf%n_v, batch_i%length, wf%n_v, batch_j%length)
            call mem%dealloc(g_iajb, batch_i%length, wf%n_v, batch_j%length, wf%n_v)
!
         enddo
      enddo
!
      wf%energy = wf%hf_energy + correlation_energy
!
   end subroutine calculate_energy_lowmem_cc2
!
!
   subroutine construct_excited_state_equation_lowmem_cc2(wf, X, R, w)
!!
!!    Construct excited state equation
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
!!    Constructs R = AX - wX, where w = X^T A X and norm(X) = sqrt(X^T X) = 1
!!
!!    Note I: we assume that X is normalized. If it is not,
!!    please normalize before calling the routine.
!!
!!    Note II: this routine constructs the excited state equation
!!    for standard CC models and the effective (!) excited state
!!    equation in perturbative models. In the lowmem_CC2 routine, for
!!    instance, X and R will be n_o*n_v vectors and A(w) will
!!    depend on the excitation energy w. See, e.g., Weigend and
!!    Hättig's RI-lowmem_CC2 paper for more on this topic. This means
!!    that w should be the previous w-value when entering the
!!    routine (so that A(w)X may be constructed approximately)
!!    in perturbative models.
!!
!!    Note III: the routine is used by the DIIS excited state solver.
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: X
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: R
!
      real(dp), intent(inout) :: w
!
      real(dp), dimension(:), allocatable :: X_copy
!
      real(dp) :: ddot
!
!     Construct residual based on previous excitation energy w
!
      call mem%alloc(X_copy, wf%n_es_amplitudes)
      call dcopy(wf%n_es_amplitudes, X, 1, X_copy, 1)
!
      call wf%effective_jacobian_transformation(w, X_copy) ! X_copy <- AX
      call dcopy(wf%n_es_amplitudes, X_copy, 1, R, 1)
      call daxpy(wf%n_es_amplitudes, -w, X, 1, R, 1)
!
!     Update excitation energy w
!
      w = ddot(wf%n_es_amplitudes, X, 1, X_copy, 1)
      call mem%dealloc(X_copy, wf%n_es_amplitudes)
!
   end subroutine construct_excited_state_equation_lowmem_cc2
!
!
end module lowmem_cc2_class
