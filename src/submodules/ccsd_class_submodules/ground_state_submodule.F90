submodule (ccsd_class) ground_state
!
!!
!!     Ground state submodule (CCSD)
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!
!!     Consists of the following module subroutines of the CCSD module:
!!
!!     ground_state_preparations:  makes preparations for the ground state solver
!!     new_amplitudes:             calculates the quasi-Newton estimate and passes the
!!                                 information needed by the DIIS routine.
!!     calc_ampeqs_norm:           calculates the norm of the amplitude equations.
!!     calc_quasi_Newton_doubles:  calculates the doubles part of the quasi-Newton estimate.
!!
!!     Can be inherited by models of the same level (e.g. CC3) without modification.
!!
!!     When inherited by higher level models (e.g. CCSDT), the new_amplitudes and calc_ampeqs_norm
!!     routines should be overridden to account for the triples quasi-Newton estimate, amplitudes,
!!     and projection vector.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine calc_ampeqs_norm_ccsd(wf, ampeqs_norm)
!!
!!    Calculate amplitude equations norm (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Calculates the norm^2 of the omega vector (i.e., the amplitude
!!    equations for the CCSD model).
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp) :: ampeqs_norm
!
      real(dp) :: ddot ! For dot product
!
      ampeqs_norm = zero
      ampeqs_norm = ddot(wf%n_t1am, wf%omega1, 1, wf%omega1, 1)
      ampeqs_norm = ddot(wf%n_t2am, wf%omega2, 1, wf%omega2, 1) + ampeqs_norm
      ampeqs_norm = sqrt(ampeqs_norm)
!
   end subroutine calc_ampeqs_norm_ccsd
!
!
   module subroutine new_amplitudes_ccsd(wf, diis_ground_state)
!
!     New amplitudes (CCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Directs the calculation of the quasi-Newton estimate Δ t_i,
!     and t_i + Δ t_i, and calls the DIIS routine to save & get
!     the amplitudes for the next iteration.
!
      implicit none
!
      class(ccsd) :: wf
!
      class(diis) :: diis_ground_state
!
      integer(i15) :: i = 0 ! for debug purposes
!
      real(dp), dimension(:,:), allocatable :: dt   ! Δ t_i
      real(dp), dimension(:,:), allocatable :: t_dt ! t_i + Δ t_i
!
!     Allocate Δ t_i and t_i + Δ t_i vectors
!
      call wf%mem%alloc(dt, wf%n_parameters, 1)
      call wf%mem%alloc(t_dt, wf%n_parameters, 1)
!
      dt   = zero
      t_dt = zero
!
!     Calculate Δ t_i
!
      call wf%calc_quasi_Newton_singles(dt)
      call wf%calc_quasi_Newton_doubles(dt)
!
!     Set t_i + Δ t_i
!
      call dcopy(wf%n_parameters, dt, 1, t_dt, 1) ! t_dt = Δ t_i
!
      call daxpy(wf%n_t1am, one, wf%t1am, 1, t_dt, 1)                   ! t_dt = t_i + Δ t_i singles
      call daxpy(wf%n_t2am, one, wf%t2am, 1, t_dt(wf%n_t1am + 1, 1), 1) ! t_dt = t_i + Δ t_i doubles
!
!     Save estimates to file and get the next amplitudes
!     (they are placed in dt on exit from diis)
!
   !   call wf%diis(dt, t_dt)
      call diis_ground_state%update(dt, t_dt, wf%disk, wf%mem)
!
!     Set the new amplitudes
!
      call dcopy(wf%n_t1am, dt, 1, wf%t1am, 1)
      call dcopy(wf%n_t2am, dt(wf%n_t1am + 1, 1), 1, wf%t2am, 1)
!
!     Deallocate vectors
!
      call wf%mem%dealloc(dt, wf%n_parameters, 1)
      call wf%mem%dealloc(t_dt, wf%n_parameters, 1)
!
   end subroutine new_amplitudes_ccsd
!
!
   module subroutine calc_quasi_Newton_doubles_ccsd(wf,dt)
!
!     Calculate quasi-Newton doubles estimate (CCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Calculates the quasi-Newton estimate Δ t_i (doubbles part)
!     and places the contribution in the dt vector (of length n_parameters,
!     with singles first, then doubles, etc. if inherited)
!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_parameters, 1) :: dt
!
      integer(i15) :: a = 0, i = 0, b = 0, j = 0
      integer(i15) :: ai = 0, bj = 0, aibj = 0, offset = 0
!
!     Calculate the doubles Δ t_i contribution
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
!                 Calculate the necessary indices
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!
                  aibj = index_packed(ai, bj)
!
                  offset = wf%n_t1am + aibj ! dt has singles first, then doubles
!
                  dt(offset,1) = - wf%omega2(aibj, 1)/(wf%fock_diagonal(wf%n_o + a, 1) + &
                                                       wf%fock_diagonal(wf%n_o + b, 1) - &
                                                       wf%fock_diagonal(i, 1) -          &
                                                       wf%fock_diagonal(j, 1))
!
               enddo
            enddo
         enddo
      enddo
!
   end subroutine calc_quasi_Newton_doubles_ccsd
!
!
   module subroutine ground_state_preparations_ccsd(wf)
!!
!!    Ground state preparations (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    A routine for preparation tasks (if any). Can be overwritten
!!    in descendants if other preparations prove necessary.
!!
      class(ccsd) :: wf
!
!     Test for the possibility of storing vir-vir-vir-vir
!     electronic repulsion integrals (g_abcd), storing the
!     integrals if possible
!
      call wf%initialize_single_amplitudes
      call wf%store_vv_vv_electronic_repulsion
!
!     Allocate double amplitudes,
!     and set the amplitudes to the MP2 guess
!
      call wf%initialize_double_amplitudes   ! Allocate double amplitudes
      call wf%construct_perturbative_doubles ! Set doubles amplitudes to MP2 guess
!
!     Allocate the projection vector
!
      call wf%initialize_omega
!
   end subroutine ground_state_preparations_ccsd
!
!
   module subroutine summary_ground_state_info_ccsd(wf, time)
!!
!!    Summary ground state info (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp) :: time
      real(dp) :: norm_singles
      real(dp) :: norm_doubles
!
      real(dp) :: t1, t2
!
!     Print energy and CPU time
!
      write(unit_output,'(/t3,a,a,a/)')'Summary of ', trim(wf%name), ' ground state calculation:'
      write(unit_output,'(t6,a25,f19.12)')  'Total energy [a.u.]:     ', wf%energy
      write(unit_output,'(t6,a25,f19.12)')  'Total time CPU (seconds): ', time
!
!     Print the dominant single excitations
!
      call print_dominant_two_index(wf%t1am, wf%n_v, wf%n_o, 'a', 'i')
!
!     Print the dominant double excitations
!
      call print_dominant_four_index(wf%t2am, wf%n_v, wf%n_o, wf%n_v, wf%n_o, &
                                                'a', 'i', 'b', 'j')
!
      write(unit_output,'(/t6,a41,f14.12/)')                &
               'Singles contribution to the full vector: ', &
               norm_singles**2/(norm_singles**2 + norm_doubles**2)
!
   end subroutine summary_ground_state_info_ccsd
!
!
end submodule ground_state
