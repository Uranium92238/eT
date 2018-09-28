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

   module subroutine omega_ccsd_a1_ccsd(wf, omega1)
!!
!!    Omega A1 term
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad, 
!!    Andreas Skeidsvoll and Alice Balbi, 2018
!!
!!    Calculates the A1 term,
!!
!!       A1: sum_ckd g_adkc * u_ki^cd,
!!
!!    and adds it to the singles projection vector (omega1) of
!!    the wavefunction object wf.
!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout):: omega1
!
!     Batching variables
!
      integer(i15) :: required        = 0
      integer(i15) :: current_a_batch = 0
!
      type(batching_index) :: batch_a
!
      real(dp), dimension(:,:), allocatable :: u_dk_ci, t_dk_ci, g_ad_kc
!
      integer(i15) :: ad_dim
!
!     u_ki^cd = 2*t_ki^cd - t_ik^cd (ordered as u_dk_ci)
!
      call mem%alloc(t_dk_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call squareup(wf%t2, t_dk_ci, (wf%n_v)*(wf%n_o))
!
      call mem%alloc(u_dk_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      u_dk_ci = - t_dk_ci
!
      call add_1432_to_1234(two, t_dk_ci, u_dk_ci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_dk_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Prepare for batching
!
!     Estimated memory required to construct g_adkc
!
      required = wf%integrals%get_required_vvov()
!
!     Initialization of the batching variable
!
      call batch_a%init(wf%n_v)                ! Initialize batching index a
      call mem%num_batch(batch_a, required) ! Determine batching information
!
!     Loop over the number of a batches
!
      do current_a_batch = 1, batch_a%num_batches
!
!        For each batch, get the limits for the a index
!
         call batch_a%determine_limits(current_a_batch)
!
         ad_dim = (batch_a%length)*(wf%n_v)
!
!        Form g_ad_kc = g_adkc
!
         call mem%alloc(g_ad_kc, ad_dim, (wf%n_o)*(wf%n_v))
!
         call wf%get_vvov(g_ad_kc,                       &
                           batch_a%first, batch_a%last,  &
                           1, wf%n_v,                    &
                           1, wf%n_o,                    &
                           1, wf%n_v)
!
         call dgemm('N','N',                     &
                     batch_a%length,             &
                     wf%n_o,                     &
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     g_ad_kc,                    & ! g_a_dkc
                     batch_a%length,             &
                     u_dk_ci,                    & ! u_dkc_i
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     omega1(batch_a%first,1),    &
                     wf%n_v)
!
         call mem%dealloc(g_ad_kc, ad_dim, (wf%n_o)*(wf%n_v))
!
      enddo ! End of batches of the index a
!
      call mem%dealloc(u_dk_ci, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine omega_ccsd_a1_ccsd
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
