submodule (ccsd_class) omega_ccsd
!
!!
!!    Omega submodule (ccsd)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Andreas Skeidsvoll and Alice Balbi, 2018
!!
!
   implicit none
!
!
!
!
contains
!
!
   module subroutine omega_ccsd_b1_ccsd(wf, omega1)
!!
!!    Omega B1
!!    Written by Sarai D. Fokestad, Eirik F. Kjønstad,
!!    Andreas Skeidsvoll and Alice Balbi, 2018
!!
!!    Calculates the B1 term,
!!
!!       B1:   - sum_ckl u_kl^ac * g_kilc,
!!
!!    and adds it to the singles projection vector (omega1) of
!!    the wavefunction object wf
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(:,:), allocatable :: g_lc_ki ! g_kilc
      real(dp), dimension(:,:), allocatable :: t_al_ck ! t_kl^ac
!     Get g_ki_lc = g_kilc
!
      call mem%alloc(g_lc_ki, (wf%n_v)*wf%n_o), (wf%n_o)**2)
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_oo(integral_type, g_lc_ki)
!
!     Form u_al_ck = u_kl^ac = 2 * t_kl^ac - t_lk^ac
!     Square up amplitudes and reorder: t_ak_cl to t_al_ck
      call mem%alloc(t_al_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*wf%n_o)
   end subroutine omega_ccsd_b1_ccsd
end submodule
