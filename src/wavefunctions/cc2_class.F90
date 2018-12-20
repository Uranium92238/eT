module cc2_class
!
!!
!!    Coupled cluster singles and perturbative doubles (cc2) 
!!    class module
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!
!
   use ccs_class
!
   implicit none
!
   type, extends(ccs) :: cc2
!
   contains
!
      procedure :: construct_omega  => construct_omega_cc2
      procedure :: omega_cc2_a1     => omega_cc2_a1_cc2
      procedure :: omega_cc2_b1     => omega_cc2_b1_cc2
      procedure :: omega_cc2_c1     => omega_cc2_c1_cc2
!
      procedure :: calculate_energy => calculate_energy_cc2
!
    !  procedure :: effective_jacobian_transformation => effective_jacobian_transformation_cc2
!
   end type cc2
!
   interface
!
      include "../submodules/cc2/omega_cc2_interface.F90"
  !    include "../submodules/cc2/jacobian_cc2_interface.F90"
!
   end interface 
!
!
contains
!
!
!
   subroutine calculate_energy_cc2(wf)
!!
!!     Calculate energy (CC2)
!!     Written by Sarai D. Folkestad, Eirik F. Kjønstad, 
!!     Andreas Skeidsvoll, 2018
!!
!!     Calculates the CC2 energy. This is only equal to the actual
!!     energy when the ground state equations are solved, of course.
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      real(dp), dimension(:,:), allocatable :: g_aibj
      real(dp), dimension(:,:), allocatable :: g_iajb
!
      real(dp) :: correlation_energy
!
      integer(i15) :: i, j, a, b, ai, bj, ia, jb, ib, ja
!
      integer(i15) :: req0, req1_i, req1_j, req2
!
      integer(i15) :: current_i_batch, current_j_batch
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
            call mem%alloc(g_aibj, wf%n_v*(batch_i%length), wf%n_v*(batch_j%length))
            call mem%alloc(g_iajb, wf%n_v*(batch_i%length), wf%n_v*(batch_j%length))
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
!$omp parallel do private(b,i,j,a,ai,ia,bj,jb,ib,ja) reduction(+:correlation_energy)
            do b = 1, wf%n_v
               do i = 1, batch_i%length
                  do j = 1, batch_j%length
                     do a = 1, wf%n_v
!
                        ai = index_two(a, i, wf%n_v)
                        ia = index_two(i, a, batch_i%length)
                        bj = index_two(b, j, wf%n_v)
                        jb = index_two(j, b, batch_j%length)
                        ib = index_two(i, b, batch_i%length)
                        ja = index_two(j, a, batch_j%length)
!
                        correlation_energy = correlation_energy + &
                                             (wf%t1(a, i + batch_i%first - 1)*wf%t1(b, j + batch_j%first - 1) &
                                             - (g_aibj(ai, bj))/(wf%fock_diagonal(wf%n_o + a, 1) &
                                                               + wf%fock_diagonal(wf%n_o + b, 1) &
                                                               - wf%fock_diagonal(i + batch_i%first - 1,1) &
                                                               - wf%fock_diagonal(j + batch_j%first - 1,1)))&
                                             *(two*g_iajb(ia,jb)-g_iajb(ib,ja))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do 
!
!
            call mem%dealloc(g_aibj, wf%n_v*(batch_i%length), wf%n_v*(batch_j%length))
            call mem%dealloc(g_iajb, wf%n_v*(batch_i%length), wf%n_v*(batch_j%length))
!
         enddo
      enddo
!
      wf%energy = wf%hf_energy + correlation_energy
!
   end subroutine calculate_energy_cc2
!
end module cc2_class
