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
   end type cc2
!
   interface
!
      include "../submodules/cc2/omega_cc2_interface.F90"
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
      real(dp) :: ddot
!
      integer(i15) :: i, j, a, b, ai, bj, ia, jb, ib, ja
!
      call mem%alloc(g_aibj, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
      call mem%alloc(g_iajb, wf%n_v*wf%n_o, wf%n_v*wf%n_o)
!
      call wf%get_vovo(g_aibj)
      call wf%get_ovov(g_iajb)
!
      wf%energy = wf%hf_energy
!
      do i = 1, wf%n_o
         do j = 1, wf%n_o
            do a = 1, wf%n_v
               do b = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
                  ia = index_two(i, a, wf%n_o)
                  bj = index_two(b, j, wf%n_v)
                  jb = index_two(j, b, wf%n_o)
                  ib = index_two(i, b, wf%n_o)
                  ja = index_two(j, a, wf%n_o)
!
                  wf%energy = wf%energy + (wf%t1(a,i)*wf%t1(b,j) &
                                       - (g_aibj(ai, bj))/(wf%fock_diagonal(wf%n_o+a,1) &
                                                         + wf%fock_diagonal(wf%n_o+b,1) &
                                                         - wf%fock_diagonal(i,1) - wf%fock_diagonal(j,1)))&
                                       *(two*g_iajb(ia,jb)-g_iajb(ib,ja))
!
               enddo
            enddo
         enddo
      enddo
!
!
   end subroutine calculate_energy_cc2
!
end module cc2_class
