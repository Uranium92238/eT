module ccsdpt_class
!
!!
!!                          CCSD(T) class module
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
!  Ancestor class module (CC3)
!
   use cc3_class
!
   implicit none 
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the CCSD(T) class -::-
!  :::::::::::::::::::::::::::::::::::::::::
!
   type, extends(cc3) :: ccsdpt
!
!     No unique variables (inherits all of them from CC3)
!
   contains
!
!     Initialization routine 
!
      procedure :: init => init_ccsdpt
!
!     Routines to construct the projection vector (omega)
!
      procedure :: construct_omega => construct_omega_ccsdpt
!
!     Energy correction routine 
!
      procedure :: calc_energy_correction => calc_energy_correction_ccsdpt
!
!     ... which is called before the ground state is destroyed:
!
      procedure :: destruct_ground_state => destruct_ground_state_ccsdpt
!
   end type ccsdpt 
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Interface to the submodule routines of CCSD(T) -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   interface 
!
!
      module subroutine construct_omega_ccsdpt(wf)
!!
!!       Construct Omega (CCSD(T))
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Directs the calculation of the projection vector (omega1, omega2)
!!       for the CCSD(T) level of theory. This is simply the calculation of 
!!       the CCSD omega vector. 
!!
         implicit none 
!
         class(ccsdpt) :: wf
!
      end subroutine construct_omega_ccsdpt
!
!
      module subroutine destruct_ground_state_ccsdpt(wf)
!!
!!       Destruct Ground State (CCSD(T))
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Calculates the CCSD(T) energy correction to the CCSD
!!       energy, thereafter deallocating the amplitudes and the 
!!       projection vector.
!!
         implicit none
!
         class(ccsdpt) :: wf
!
      end subroutine destruct_ground_state_ccsdpt
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
   subroutine init_ccsdpt(wf)
!!
!!    Initialize CCSD(T) object
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
      class(ccsdpt) :: wf
!
!     Set model name 
!
      wf%name = 'CCSD(T)'
!
!     Set implemented methods
!
      wf%implemented%ground_state = .true.
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
   end subroutine init_ccsdpt
!
!
   subroutine calc_energy_correction_ccsdpt(wf)
!!
!!    Calculate Energy Correction (CCSD(T))
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Calculates the CCSD(T) correction to the CCSD energy:
!!
!!       E(CCSD(T)) = E(CCSD) + sum_ai u_ai v_ai + sum_aibj u_ai_bj v_ai_bj,
!!
!!    where
!!
!!       v_ai    = sum_cdkl (t_ikl^acd - t_lki^acd) L_kcld 
!!       v_ai_bj = sum_cdk (t_ijk^acd L_bckd - t_kji^acd g_kdbc) -
!!                 sum_ckl (t_ikl^abc L_kjlc - t_lki^abc g_kjlc)
!!
!!    and
!!
!!       u_ai    = 2 * t_i^a 
!!       u_ai_bj = 4 t_ij^ab - 2 t_ji^ab. 
!!
!!    Note: the vector v is similar to the CC3 omega vector. The 
!!    routine therefore relies heavily on the omega routines of 
!!    the parent CC3 class.
!!
!!    Note: as is the case for CC3 omega, this routine (in particular
!!    the i >= j >= k saving) is not optimized. Both CC3 and CCSD(T)
!!    should be optimized simultaneously.
!!
      class(ccsdpt) :: wf 
!
      real(dp), dimension(:,:), allocatable :: u_ai
      real(dp), dimension(:,:), allocatable :: u_ai_bj 
!
      real(dp), dimension(:,:), allocatable :: v_ai
      real(dp), dimension(:,:), allocatable :: v_ai_bj
!
      integer(i15) :: a = 0, i = 0, j = 0, k = 0, b = 0, ai = 0 
      integer(i15) :: bj = 0, aj = 0, bi = 0, aibj = 0, ajbi = 0
!
      real(dp), dimension(:,:), allocatable :: t_abc ! Triples, t_ijk^abc, fixed ijk
!
      real(dp) :: ddot, correction 
!
!     Let the user know the CCSD(T) correction is being computed
!
      write(unit_output,'(t3,a/,t3,a/)') 'Adding the CCSD(T) correction to the above,',& 
                                         'iterated CCSD energy.'
      flush(unit_output)
!
!     :: Calculate the v vector ::
!
      call wf%mem%alloc(v_ai, wf%n_v, wf%n_o)
!
!     Write to disk the integrals needed for calculation of v
!
!     Note: These are the same as for CC3 omega, without the
!           T1 transformation. We therefore reuse this routine.
!
      v_ai = wf%t1am          ! Make a copy of the singles
      wf%t1am = zero          ! Set the singles to zero 
!
      call wf%omega_integrals ! Get the (non-T1) transformed integrals
! 
      wf%t1am = v_ai          ! Put the singles amplitudes back
      v_ai = zero             ! Reset the W_ai vector 
!
      wf%omega1 = zero        ! Will hold the singles part of W, 
                              ! arising from the E1 term 
!
      call wf%mem%alloc(t_abc, (wf%n_v)**3, 1)
      t_abc = zero
!
      call wf%mem%alloc(v_ai_bj, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      v_ai_bj = zero 
!
      do i = 1, wf%n_o
         do j = 1, wf%n_o
            do k = 1, wf%n_o
!
!              Calculate the triples amplitudes:
!
!              Calculate W_abc = P_ijk^abc ( sum_d t_ij^ad g_bdck - sum_l t_il^ab g_ljck )
!              and divide by orbital energy difference, t_abc = - W_abc / e_abc
!
               t_abc = zero
               call wf%calc_triples(t_abc,i,j,k)         
!
!              Note: v_ai is identical to the CC3 omega1(a,i), 
!              with non-transformed integrals:
!
               call wf%omega_cc3_a1(t_abc,i,j,k)
!
!              Note: v_ai_bj is identical the CC3 B2 contribution
!              to omega2(ai,bj):
!
               call wf%omega_cc3_b2(v_ai_bj,t_abc,i,j,k)
!
            enddo
         enddo
      enddo
!
      call dcopy((wf%n_o)*(wf%n_v), wf%omega1, 1, v_ai, 1)
!
!     Deallocate the triples amplitude 
!
      call wf%mem%dealloc(t_abc, (wf%n_v)**3, 1)
!
!     Calculate the u vector 
!
      call wf%mem%alloc(u_ai, wf%n_v, wf%n_o)
      u_ai = zero 
!
      call wf%mem%alloc(u_ai_bj, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      u_ai_bj = zero 
!
      call daxpy((wf%n_o)*(wf%n_v), two, wf%t1am, 1, u_ai, 1) ! u_ai = 2 * t_ai 
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = index_two(a, i, wf%n_v)
!
            do j = 1, wf%n_o
!
               aj = index_two(a, j, wf%n_v)
!
               do b = 1, wf%n_v
!
                  bj = index_two(b, j, wf%n_v)
                  bi = index_two(b, i, wf%n_v)
!
                  aibj = index_packed(ai, bj)
                  ajbi = index_packed(aj, bi)
!
                  u_ai_bj(ai,bj) = four*(wf%t2am(aibj,1)) - two*(wf%t2am(ajbi,1))
!
               enddo
            enddo
         enddo
      enddo
!
!     Find the energy corrections and print the CCSD(T) energy
!     to the main output file
!
      correction = zero
      correction = ddot((wf%n_o)*(wf%n_v), u_ai, 1, v_ai, 1)
!
      wf%energy = wf%energy + correction  
!
      correction = zero 
      correction = ddot((wf%n_o)*(wf%n_v)*(wf%n_o)*(wf%n_v), u_ai_bj, 1, v_ai_bj, 1)
!
      wf%energy = wf%energy + correction
!
      write(unit_output,'(t3,a27,f14.8)') 'CCSD(T) energy (hartrees):', wf%energy
!
!     Deallocations 
!
      call wf%mem%dealloc(u_ai, wf%n_v, wf%n_o)
      call wf%mem%dealloc(u_ai_bj, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call wf%mem%dealloc(v_ai, wf%n_v, wf%n_o)
      call wf%mem%dealloc(v_ai_bj, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      end subroutine calc_energy_correction_ccsdpt
!
!
end module ccsdpt_class