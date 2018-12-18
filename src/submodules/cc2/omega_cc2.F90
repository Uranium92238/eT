submodule (cc2_class) omega_cc2
!
!!
!!    Omega submodule (cc2)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Linda Goletto, 2018
!!
!!    Routines to construct 
!!
!!    Ω =  < mu | exp(-T) H exp(T) | R >
!!
!
   implicit none
!
!
contains
!
   module subroutine construct_omega_cc2(wf, omega)
!!
!!    Construct omega 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Linda Goletto, 2018
!!
!!    Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current amplitudes of the object wfn
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
      omega = zero
!
      call wf%omega_ccs_a1(omega)
!
      call wf%omega_cc2_a1(omega)
      call wf%omega_cc2_b1(omega)
      call wf%omega_cc2_c1(omega)
!
   end subroutine construct_omega_cc2
!
!
   module subroutine omega_cc2_a1_cc2(wf, omega)
!!
!!    Omega CC2 A1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Linda Goletto, 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
   end subroutine omega_cc2_a1_cc2
!
!
   module subroutine omega_cc2_b1_cc2(wf, omega)
!!
!!    Omega CC2 B1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Linda Goletto, 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
   end subroutine omega_cc2_b1_cc2
!
!
   module subroutine omega_cc2_c1_cc2(wf, omega)
!!
!!    Omega CC2 C1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Linda Goletto, 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
      real(dp), dimension(:,:), allocatable :: g_aibj
      real(dp), dimension(:,:), allocatable :: L_aijb
!
      integer(i15) :: i, j, a, b
      integer(i15) :: ai, aj, bi, bj, jb
!
      call mem%alloc(g_aibj, wf%n_o*wf%n_v, wf%n_o*wf%n_v)
!
      call wf%get_vovo(g_aibj)
!
      call mem%alloc(L_aijb, wf%n_o*wf%n_v, wf%n_o*wf%n_v)
!
      do i = 1, wf%n_o 
         do j = 1, wf%n_o
            do a = 1, wf%n_v
!
               ai = wf%n_v*(i-1) + a
               aj = wf%n_v*(j-1) + a
!
               do b = 1, wf%n_v

!
                  bi = wf%n_v*(i-1) + b
                  bj = wf%n_v*(j-1) + b
                  jb = wf%n_o*(b-1) + j
!
                  L_aijb(ai, jb) = -two*(g_aibj(ai, bj)/(wf%fock_ab(a,a) + wf%fock_ab(b,b) &
                                                      - wf%fock_ij(i,i) - wf%fock_ij(j,j))) &
                                    + (g_aibj(aj, bi)/(wf%fock_ab(a,a) + wf%fock_ab(b,b) &
                                                      - wf%fock_ij(i,i) - wf%fock_ij(j,j))) 
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_aibj, wf%n_o*wf%n_v, wf%n_o*wf%n_v)
!
      call dgemm('N', 'N', &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  one,           &
                  L_aijb,        &
                  wf%n_v*wf%n_o, &
                  wf%fock_ia,    &
                  wf%n_v*wf%n_o, &
                  one,           &
                  omega,         &
                  wf%n_v*wf%n_o)
!
      call mem%dealloc(L_aijb, wf%n_o*wf%n_v, wf%n_o*wf%n_v)
!
   end subroutine omega_cc2_c1_cc2
!
end submodule