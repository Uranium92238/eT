module mp2_class
!
!!
!!              Second-order Möller-Plesset (MP2) class module                                 
!!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017         
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
!  The ancestor class module (HF)
!
   use hf_class
!
   implicit none 
!
!  :::::::::::::::::::::::::::::::::::::
!  -::- Definition of the MP2 class -::-
!  ::::::::::::::::::::::::::::::::::::: 
!
   type, extends(hf) :: mp2
!
!     No unique variables (inherits all of them from HF)
!
   contains
!
!     Initialization and driver routines
!
      procedure :: init => init_mp2
      procedure :: drv  => drv_mp2
!
!     Routine to calculate the MP2 energy
!
      procedure :: calc_energy => calc_energy_mp2
!
   end type mp2
!
!
contains
!
!  ::::::::::::::::::::::::::::::::::::::::::::
!  -::- Initialization and driver routines -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::
!
   subroutine init_mp2(wf)
!!
!!    Initialize MP2 object
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Performs the following tasks
!!
!!    - Sets HF orbital and energy information by reading from file
!!    - Transforms AO Cholesky vectors to MO basis and saves to file 
!!
      implicit none 
!
      class(mp2) :: wf
!
!     Set model name 
!
      wf%name = 'MP2'
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
   end subroutine init_mp2
!
!
   subroutine drv_mp2(wf)
!!
!!    MP2 Driver
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    The driver for MP2 checks whether the ground state is requested,
!!    and if so, calls the calculation of the MP2 energy.
!!
      implicit none 
!
      class(mp2) :: wf
!
      if (wf%tasks%ground_state) then 
!
!        Ground state calculation requested
!
         if (wf%implemented%ground_state) then 
!
            call wf%calc_energy
!
         else
!
            write(unit_output,'(t3,a,a)') &
               'Error: ground state solver not implemented for ',trim(wf%name)
            stop
!
         endif
!
      else
!
         write(unit_output,'(t3,a,a)') &
            'Error: only the ground state is implemented for ',trim(wf%name)
         stop
!
      endif
!
   end subroutine drv_mp2
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Class subroutines and functions -::- 
!  :::::::::::::::::::::::::::::::::::::::::
!
   subroutine calc_energy_mp2(wf)
!!
!!    Calculate Energy (MP2)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017 
!!
!!    Calculates the MP2 energy. 
!!
      implicit none 
!
      class(mp2) :: wf 
!
      real(dp), dimension(:,:), allocatable :: L_ia_J  ! L_ia^J
      real(dp), dimension(:,:), allocatable :: g_ia_jb ! g_iajb
!
      integer(i15) :: i = 0, a = 0, j = 0, b = 0
      integer(i15) :: ia = 0, jb = 0, ib = 0, ja = 0
!
      real(dp) :: begin_timer, end_timer, integral_time
!
      write(unit_output,'(/t3,a/)') 'Calculating the MP2 energy:'
      flush(unit_output)
!
      call cpu_time(begin_timer)
!
!     Read L_ia^J from disk 
!
      call wf%mem%alloc(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call wf%read_cholesky_ia(L_ia_J)
!
!     Form g_iajb = sum_J L_ia^J L_ia^J
!
      call wf%mem%alloc(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ia_jb,           &
                  (wf%n_o)*(wf%n_v))
!
      call wf%mem%dealloc(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Add the MP2 correction to the Hartree-Fock energy 
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ia = index_two(i, a, wf%n_o)
!
            do b = 1, wf%n_v
!
               ib = index_two(i, b, wf%n_o)
!
               do j = 1, wf%n_o
!
                  ja = index_two(j, a, wf%n_o)
                  jb = index_two(j, b, wf%n_o)
!
                  wf%energy = wf%energy - &
                           (g_ia_jb(ia,jb)*(two*g_ia_jb(ia,jb)-g_ia_jb(ib,ja)))/   &
                                                   (wf%fock_diagonal(wf%n_o+a,1) + &
                                                    wf%fock_diagonal(wf%n_o+b,1) - &
                                                    wf%fock_diagonal(i,1) -        &
                                                    wf%fock_diagonal(j,1))
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Print the energy to the main output file
!
      write(unit_output,'(t3,a27,f14.8)') 'Total energy (hartrees):', wf%energy
!
      call cpu_time(end_timer)
      write(unit_output,'(t3,a27,f14.8/)') 'Total time (seconds):', end_timer-begin_timer
      flush(unit_output)
!
   end subroutine calc_energy_mp2
!
!
end module mp2_class