module uhf_class
!
!!
!!    Unrestricted Hartree-Fock (HF) class module
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!
   use hf_class
!
   use reordering
   use array_utilities
   use array_analysis
   use interval_class
   use index
!
   implicit none
!
   type, extends(hf):: uhf
!
      integer(i15) :: n_alpha
      integer(i15) :: n_beta 
!
      real(dp) :: uhf_energy
!
      real(dp), dimension(:,:), allocatable :: ao_density_a 
      real(dp), dimension(:,:), allocatable :: ao_density_b
! 
      real(dp), dimension(:,:), allocatable :: ao_fock_a 
      real(dp), dimension(:,:), allocatable :: ao_fock_b
!
      real(dp), dimension(:,:), allocatable :: orbital_coefficients_a
      real(dp), dimension(:,:), allocatable :: orbital_coefficients_b      
!
      real(dp), dimension(:,:), allocatable :: orbital_energies_a
      real(dp), dimension(:,:), allocatable :: orbital_energies_b
!
	contains
!
      procedure :: prepare => prepare_uhf
!
      procedure :: determine_n_alpha_and_beta_electrons => determine_n_alpha_and_beta_electrons_uhf
!
      procedure :: construct_ao_spin_density => construct_ao_spin_density_uhf
      procedure :: construct_ao_spin_fock    => construct_ao_spin_fock_uhf
!
      procedure :: form_ao_density => form_ao_density_uhf
!
   end type uhf
!
!
contains 
!
!
   subroutine prepare_uhf(wf)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      wf%name = 'UHF'
!
      call wf%system%prepare()
!
      call get_n_aos(wf%n_ao)
!
      call initialize_coulomb()
      call initialize_kinetic()
      call initialize_nuclear()
      call initialize_overlap()
!
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
      call wf%decompose_ao_overlap() 
!
      wf%n_o = (wf%system%get_n_electrons())/2
      wf%n_v = wf%n_mo - wf%n_o
!
      call wf%determine_n_alpha_and_beta_electrons()
!
   end subroutine prepare_uhf
!
!
   subroutine determine_n_alpha_and_beta_electrons_uhf(wf)
!!
!!    Determine the number of alpha and beta electrons 
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
!!    Determines the number of alpha and beta elctrons 
!!    from the provided multiplicity and the number of
!!    electrons in the system. We assume that n_alpha 
!!    is greater than or equal to n_beta without loss
!!    of generality.
!!
      implicit none 
!
      class(uhf) :: wf 
!
      integer(i15) :: n_alpha_m_n_beta, n_alpha_p_n_beta 
!
      n_alpha_m_n_beta = 2*(wf%system%multiplicity)
      n_alpha_p_n_beta = wf%system%n_electrons 
!
      wf%n_alpha = (n_alpha_p_n_beta + n_alpha_m_n_beta)/2
      wf%n_beta  = wf%system%n_electrons - wf%n_alpha
!
      if ((wf%n_alpha + wf%n_beta) .ne. wf%system%n_electrons) then 
!
         call output%error_msg('Given multiplicity and number of electrons is inconsistent.')
!
      endif 
!
   end subroutine determine_n_alpha_and_beta_electrons_uhf
!
!
   subroutine construct_ao_spin_density_uhf(wf, sigma)
!!
!!    Construct AO spin density 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Constructs either the alpha or beta density, 
!!    depending on the value of 'spin' (='alpha' or 'beta').
!!    The routine assumes that the associated orbital coefficients 
!!    are allocated and set. Typically, the construction of these
!!    densities will follow a Roothan-Hall step. 
!!
!!    We employ the normalization so that D = D_alpha + D_beta is 
!!    equal to the restricted HF density if the alpha and beta 
!!    densities are equal.
!!
      implicit none 
!
      class(uhf) :: wf 
!
      character(len=*), intent(in) :: sigma
!
      if (trim(sigma) == 'alpha') then 
!
         call dgemm('N', 'T',                   &
                     wf%n_ao,                   &
                     wf%n_ao,                   &
                     wf%n_o,                    &
                     one,                       &
                     wf%orbital_coefficients_a, &
                     wf%n_ao,                   &
                     wf%orbital_coefficients_a, &
                     wf%n_ao,                   &
                     zero,                      &
                     wf%ao_density_a,           &
                     wf%n_ao)
!
      elseif (trim(sigma) == 'beta') then 
!
         call dgemm('N', 'T',                   &
                     wf%n_ao,                   &
                     wf%n_ao,                   &
                     wf%n_o,                    &
                     one,                       &
                     wf%orbital_coefficients_b, &
                     wf%n_ao,                   &
                     wf%orbital_coefficients_b, &
                     wf%n_ao,                   &
                     zero,                      &
                     wf%ao_density_b,           &
                     wf%n_ao)
!
      else
!
         call output%error_msg('Did not recognize spin variable in construct_ao_spin_density:' // trim(sigma))
!
      endif  
!
   end subroutine construct_ao_spin_density_uhf
!
!
   subroutine form_ao_density(wf)
!!
!!    Form AO density 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Given that the spin densities have been constructed,
!!    this routine sets 
!!
!!       D = D_alpha + D_beta     
!!
!!    Note that it is not constructed from the orbital
!!    coefficients and should be preceded by calls to 
!!    "construct_ao_spin_density" to make sure that the 
!!    spin densities are up-to-date.
!!
      implicit none 
!
      class(uhf) :: wf 
!
      wf%ao_density = wf%ao_density_a + wf%ao_density_b
!
   end subroutine form_ao_density
!
!
   subroutine construct_ao_spin_fock_uhf(wf, D, D_sigma, sigma)
!!
!!    Construct AO spin Fock 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    The routine computes the alpha or beta Fock matrix, depending 
!!    on the value of the spin 'sigma' (='alpha' or 'beta').  
!!
!!    Note that this routine has no screening as of today,
!!    given that UHF is mainly useful in SAD (for us, today), 
!!    where efficiency really isn't a big concern. It uses 
!!    full permutational symmetry, like its more efficient 
!!    variant in the restricted HF class.
!!
      implicit none 
!
      class(uhf) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D 
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D_sigma 
!
      character(len=*), intent(in) :: sigma 
!
      if (trim(sigma) == 'alpha') then 
!
         
!
      elseif (trim(sigma) == 'beta') then 
!
      else 
!
         call output%error_msg('Did not recognize spin variable in construct_ao_fock:' // trim(sigma))
!
      endif 
!
   end subroutine construct_ao_spin_fock_uhf
!
!
end module uhf_class
