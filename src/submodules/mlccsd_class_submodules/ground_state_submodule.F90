submodule (mlccsd_class) ground_state
!
!!    
!!     Ground state submodule (MLCCSD)
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!    
!!    
!!     Consists of the following module subroutines of the MLCCSD module:
!!     
!!     new_amplitudes:             Calculates the quasi-Newton estimate and passes the 
!!                                 information needed by the DIIS routine.
!!     calc_ampeqs_norm:           Calculates the norm of the amplitude equations.
!!     calc_quasi_Newton_doubles:  Calculates the doubles part of the quasi-Newton estimate.
!!     initialize_ground_state:    Initializes the amplitudes (MP2 estimate) and the amplitude 
!!                                 equations.
!!    
!!    
!
   implicit none 
!
!
contains
!
!
   module subroutine calc_ampeqs_norm_mlccsd(wf, ampeqs_norm)
!
!     Calculate Amplitude Equations Norm (MLCCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
      implicit none 
!
      class(mlccsd) :: wf 
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
   end subroutine calc_ampeqs_norm_mlccsd
!
!
   module subroutine new_amplitudes_mlccsd(wf)
!
!     New Amplitudes (MLCCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Directs the calculation of the quasi-Newton estimate Δ t_i, 
!     and t_i + Δ t_i, and calls the DIIS routine to save & get 
!     the amplitudes for the next iteration.
!
      implicit none 
!
      class(mlccsd) :: wf 
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
      call wf%diis(dt, t_dt)
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
   end subroutine new_amplitudes_mlccsd
!
!
   module subroutine calc_quasi_Newton_doubles_mlccsd(wf,dt)
!
!     Calculate quasi-Newtoni doubles estimate (CCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Calculates the quasi-Newton estimate Δ t_i (doubbles part)
!     and places the contribution in the dt vector (of length n_parameters,
!     with singles first, then doubles, etc. if inherited)
!
      implicit none 
!
      class(mlccsd) :: wf 
!
      real(dp), dimension(wf%n_parameters, 1) :: dt
!
      integer(i15) :: a = 0, i = 0, b = 0, j = 0
      integer(i15) :: ai = 0, bj = 0, aibj = 0, offset = 0
!
!     Active space variables
!
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
!
!     Calculate first/last indeces
! 
      call wf%get_CCSD_active_indices(first_active_o, first_active_v)
      call wf%get_CCSD_n_active(n_active_o, n_active_v)
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1 
!
!     Calculate the doubles Δ t_i contribution
!
      do a = 1, n_active_v
         do i = 1, n_active_o
            do b = 1, n_active_v
               do j = 1, n_active_o
!
!                 Calculate the necessary indices 
!
                  ai = index_two(a, i, n_active_v)
                  bj = index_two(b, j, n_active_v)
!
                  aibj = index_packed(ai, bj) 
!
                  dt(wf%n_t1am + aibj ,1) = - wf%omega2(aibj, 1)/&
                                    (wf%fock_diagonal(wf%n_o + a + first_active_v - 1, 1) + &
                                     wf%fock_diagonal(wf%n_o + b + first_active_v - 1, 1) - &
                                     wf%fock_diagonal(i + first_active_o - 1, 1) -          &
                                     wf%fock_diagonal(j + first_active_o - 1, 1))
!
               enddo
            enddo
         enddo
      enddo

!
   end subroutine calc_quasi_Newton_doubles_mlccsd
!
!
   module subroutine initialize_ground_state_mlccsd(wf)
!!
!!    Initialize Ground State (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Initializes the amplitudes and the projection vector for the ground
!!    state solver.
!!
      implicit none 
!
      class(mlccsd) :: wf
!
      call wf%initialize_amplitudes          ! Allocate amplitudes
      call wf%construct_perturbative_doubles ! Set doubles amplitudes to MP2 guess 
      call wf%initialize_omega               ! Allocate projection vector 
!
   end subroutine initialize_ground_state_mlccsd
!
!
end submodule ground_state 