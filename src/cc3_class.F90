module cc3_class
!
!!
!!                           CC3 class module
!!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
!  General tools
!
   use types
   use utils
   use workspace
   use input_output
!
!  Ancestor class module (CCSD)
!
   use ccsd_class
!
   implicit none 
!
!  :::::::::::::::::::::::::::::::::::::
!  -::- Definition of the CC3 class -::-
!  :::::::::::::::::::::::::::::::::::::
!
   type, extends(ccsd) :: cc3
!
!     No unique variables (inherits all of them from CCSD)
!
   contains
!
!     Initialization routine 
!
      procedure :: init => init_cc3
!
!     Routines to construct the projection vector (omega)
!
      procedure :: construct_omega => construct_omega_cc3
!
!     Helper routines 
!
      procedure :: omega_integrals => omega_integrals_cc3
      procedure :: calc_triples    => calc_triples_cc3
!
      procedure :: omega_cc3_a1 => omega_cc3_a1_cc3
!
      procedure :: omega_cc3_a2 => omega_cc3_a2_cc3
      procedure :: omega_cc3_b2 => omega_cc3_b2_cc3
!
   end type cc3 
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Interface to the submodule routines of CC3 -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   interface 
!
!
      module subroutine construct_omega_cc3(wf)
!!
!!       Construct Omega (CC3)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Directs the calculation of the projection vector (omega1, omega2)
!!       for the CC3 level of theory.
!!
         implicit none 
!
         class(cc3) :: wf
!
      end subroutine construct_omega_cc3
!
!
      module subroutine omega_integrals_cc3(wf)
!!
!!       Omega Integrals (CC3)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Calculates g_bdck and saves to disk in the order bcd, k (one record per d and k)
!!       Calculates g_ljck and saves to disk in the order lc, jk (one record per j and k)
!!       Calculates g_jbkc and saves to disk in the order bc, kj (one record per k and j)
!!       Calculates g_ilkc and saves to disk in the order cl, ik (one record per i and k)
!!       Calculates g_dbkc and saves to disk in the order bcd, k (one record per d and k)
!!
         implicit none 
!
         class(cc3) :: wf 
!
      end subroutine omega_integrals_cc3
!
!
      module subroutine calc_triples_cc3(wf,w_abc,i,j,k)
!!
!!       Calculate Triples (CC3)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Calculate W_abc = P_ijk^abc ( sum_d t_ij^ad g_bdck - sum_l t_il^ab g_ljck )
!!       and divide by orbital energy difference, t_abc = - W_abc / e_abc. On exit,
!!       w_abc = t_abc, the CC3 triples amplitude t_ijk^abc for a given set of occupied
!!       indices i, j, and k.
!!
         implicit none 
!
         class(cc3) :: wf 
!
         real(dp), dimension((wf%n_v)**3, 1) :: w_abc
!
         integer(i15), intent(in) :: i, j, k
!
      end subroutine calc_triples_cc3
!
!
      module subroutine omega_cc3_a1_cc3(wf,t_abc,i,j,k)
!!
!!       Omega A1 (CC3)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Calculates the A1 term,
!! 
!!          sum_bc (t_ijk^abc - t_ijk^cba) L_jbkc,
!!
!!       for a given set of i, j, and k, and adds the contribution to 
!!       the singles projection vector (omega1).
!!
         implicit none 
!
         class(cc3) :: wf 
!
         real(dp), dimension((wf%n_v)**3, 1) :: t_abc 
!
         integer(i15), intent(in) :: i, j, k
!
      end subroutine omega_cc3_a1_cc3
!
!
      module subroutine omega_cc3_a2_cc3(wf,omega_ai_bj,t_abc,i,j,k)
!!
!!       Omega A2 (CC3)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Calculates the A2 term,
!! 
!!          sum_c (t_ijk^abc - t_ijk^cba) F_kc,
!!
!!       for a given set of i, j, and k, and adds the contribution to 
!!       the doubles projection vector (omega2), element ai_bj.
!!
         implicit none 
!
         class(cc3) :: wf 
!
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: omega_ai_bj 
         real(dp), dimension((wf%n_v)**3, 1) :: t_abc 
!
         integer(i15), intent(in) :: i, j, k
!
      end subroutine omega_cc3_a2_cc3
!
!
      module subroutine omega_cc3_b2_cc3(wf,omega_ai_bj,t_abc,i,j,k)
!!
!!       Omega B2 (CC3)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Calculates the B2 term,
!! 
!!          omega(al,bj) = - sum_c  (2 t_ijk^abc - t_ijk^cba - t_ijk^acb) g_ilkc
!!          omega(ai,dj) = + sum_bc (2 t_ijk^abc - t_ijk^cba - t_ijk^acb) g_dbkc
!!
!!       for a given set of i, j, and k, and adds the contribution to 
!!       the doubles projection vector (omega2).
!!
         implicit none 
!
         class(cc3) :: wf 
!
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: omega_ai_bj 
         real(dp), dimension((wf%n_v)**3, 1) :: t_abc 
!
         integer(i15), intent(in) :: i, j, k
!
      end subroutine omega_cc3_b2_cc3
!
!
   end interface
!
!
contains
!
!
!  ::::::::::::::::::::::::::::::::::::::::::::
!  -::-       Initialization routine       -::-
!  ::::::::::::::::::::::::::::::::::::::::::::
!
!
   subroutine init_cc3(wf)
!!
!!    Initialize CC3 object
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Performs the following tasks:
!!
!!       1. Sets HF orbital and energy information by reading from file (read_hf_info)
!!       2. Transforms AO Cholesky vectors to MO basis and saves to file (read_transform_cholesky)
!!       3. Allocates the Fock matrix and sets it to zero
!!       4. Initializes the amplitudes (sets their initial values and associated variables)
!!
!!    Note: this routine does not calculate the energy, which is postponed until the wavefunction
!!    is passed to the ground-state solver.
!!
      implicit none 
!
      class(cc3) :: wf
!
      integer(i15) :: unit_input = -1
!
!     Set model name 
!
      wf%name = 'CC3'

!     Open input file eT.inp
!
      call generate_unit_identifier(unit_input)
      open(unit=unit_input, file='eT.inp', status='old', form='formatted')
      rewind(unit_input)
!
!     Read general specifications (memory and diskspace for calculation)
!
      call wf%general_specs_reader(unit_input)
!
!     Set implemented methods
!
      wf%implemented%ground_state = .true.
!
!     Read calculation tasks from input file eT.inp
!     
      call wf%calculation_reader(unit_input)
!
!     Close input file
!
      close(unit_input)
!
!     Read Hartree-Fock info from SIRIUS
!
      call wf%read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf%read_transform_cholesky 
!
!     Initialize (singles and doubles) amplitudes
!
      call wf%initialize_amplitudes
!
!     Set the number of parameters in the wavefunction
!     (that are solved for in the ground and excited state solvers) 
!
      wf%n_parameters = wf%n_t1am + wf%n_t2am
!
!     Initialize the Fock matrix (allocate and construct given the initial amplitudes)
!
      call wf%initialize_fock_matrix
!
   end subroutine init_cc3
!
!
end module cc3_class
