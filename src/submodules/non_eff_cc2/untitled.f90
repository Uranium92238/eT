submodule (non_eff_cc2_class) omega_cc2
!
!!
!!    Omega submodule (CC2)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto and Alexander Paul, 2018
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
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Direqts the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current wavefunction amplitudes.
!!
      implicit none
!
      class(non_eff_cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
      real(dp), dimension(:,:,:,:), allocatable :: u_aibj
!
      omega = zero
!
      call wf%omega_ccs_a1(omega)
!
      call mem%alloc(u_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%construct_u(u_aibj)
!
      call wf%omega_cc2_a1(omega, u_aibj, omega)
      call wf%omega_cc2_b1(omega, u_aibj, omega)
      call wf%omega_cc2_c1(omega, u_aibj, omega)
!
      call mem%dealloc(u_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine construct_omega_cc2
!
!
   module subroutine omega_cc2_a1_cc2(wf, u_bjci, omega)
!!
!!    Omega CC2 A1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Calculates the A1 term,
!!
!!       A1: sum_ckd u_bj_ci * g_abjc,
!!
!!    with 
!!       
!!       u_bj_ci = 2*t_bj_ci - t_bi_cj
!!
!!    and adds it to the projection vector omega
!!
      implicit none
!
      class(non_eff_cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u_bjci
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout)          :: omega
!
   end subroutine omega_cc2_a1_cc2
!
!
   module subroutine omega_cc2_b1_cc2(wf, u_ajbk, omega)
!!
!!    Omega CC2 B1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Calculates the B1 term,
!!
!!       B1: - sum_ckl g_kb,ji * u_aj,bk,
!!
!!    with
!!
!!      u_aj_bk = 2t_aj,bk - t_ak,bj
!!
!!    and adds it to the projection vector (omega) of
!!    the wavefunction object wf.
!!
      implicit none
!
      class(non_eff_cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u_ajbk
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout)          :: omega
!
   end subroutine omega_cc2_b1_cc2
!
!
   module subroutine omega_cc2_c1_cc2(wf, u_aibj, omega)
!!
!!    Omega CC2 C1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,  Jan 2019
!!
!!    Calculates the C1 term,
!!
!!       C1: - sum_bj u_ai,bj * F_jb,
!!
!!    with 
!!       
!!       u_ai_bj = 2*t_ai_bj - t_aj_bi
!!
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u_aibj
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout)          :: omega
!
    end subroutine omega_cc2_c1_cc2
!
end submodule omega_cc2
