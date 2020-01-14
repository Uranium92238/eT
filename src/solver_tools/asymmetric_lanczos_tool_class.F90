!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
!
module asymmetric_lanczos_tool_class
!
!!
!!    Asymmetric Lanczos solver tool class
!!    Written by Torsha Moitra and Sarai D. Folkestad 
!!    and Sonia Coriani, Sep-Nov 2019
!!
!!    Constructs the tridiagonal matrix
!!
!!       T = P^T A Q
!!
!!    which is the asymmetric transformation 
!!    of some matrix A. P and Q satisfies
!!    P^T Q = I.
!!    
!!    T is diagonalized to yield left and right
!!    eigenvectors in addition to the eigenvalues.
!!
!!    We have
!!
!!       L^T T R = D = L^T P^T A Q R
!!
!!    The diagonal elements of T are 
!!    given as 
!!
!!       T(i,i) = alpha_i = p_i^T A q_i
!!
!!    The elements below/above the diagonal 
!!    are defined as
!!
!!       T(i,i-1) = beta_i-1 = sqrt(abs(s_i-1^T r_i-1))
!!       T(i-1,i) = gamma_i-1 = sgn(s_i-1^T r_i-1) beta_i-1
!!
!!    where
!!
!!       s_i = p_i^T A - beta_i-1 p_i-1 - alpha_i p_i
!!       r_i = A q_i - gamma_i-1 q_i-1 - alpha_i q_i.
!!
!!    We have 
!!
!!       p_i+1 = s_i/gamma_i
!!       q_i+1 = r_i/beta_i
!!
!!
!!    See for further details: 
!!       S. Coriani et. al., J. Chem. Theory Comput. 2012, 8, 5, 1616-1628
!!  
!  
use kinds
use parameters
use sequential_file_class, only: sequential_file
use memory_manager_class, only: mem
use global_out, only: output
!
   type :: asymmetric_lanczos_tool
!
      integer :: chain_length
      integer :: n_parameters
!
      character (len=200) :: normalization
!
      real(dp), dimension (:), allocatable :: alpha_, beta_, gamma_ 
!
      type(sequential_file), allocatable :: p_file, q_file, T_file
!
   contains
!
      procedure :: read_p                    => read_p_asymmetric_lanczos_tool
      procedure :: read_q                    => read_q_asymmetric_lanczos_tool
!
      procedure :: calculate_alpha           => calculate_alpha_asymmetric_lanczos_tool
      procedure :: construct_s               => construct_s_asymmetric_lanczos_tool
      procedure :: construct_r               => construct_r_asymmetric_lanczos_tool
      procedure :: calculate_beta_gamma_p_q  => calculate_beta_gamma_p_q_asymmetric_lanczos_tool
      procedure :: diagonalize_T             => diagonalize_T_asymmetric_lanczos_tool
      procedure :: save_T                    => save_T_asymmetric_lanczos_tool
      procedure :: biorthonormalize          => biorthonormalize_asymmetric_lanczos_tool
      procedure :: cleanup                   => cleanup_asymmetric_lanczos_tool
!
   end type asymmetric_lanczos_tool
!
   interface asymmetric_lanczos_tool
!
      procedure :: new_asymmetric_lanczos_tool
!
   end interface asymmetric_lanczos_tool

contains
!
!
   function new_asymmetric_lanczos_tool(n_parameters, chain_length, p1, q1, normalization) &
                    result(lanczos)
!!
!!    New asymmetric Lanczos tool
!!    Written by Torsha Moitra, Sarai D. Folkestad 
!!    and Sonia Coriani, Sep-Nov 2019
!!
!!    Constructor for the asymmetric Lanczos tool.
!!
!!    Allocates arrays (alpha_, beta_, gamma_)
!!    and initializes files.
!!
!!    "chain_length" : dimension of the Krylov subspace
!!    
!!    "p_1" : initial left vector
!!    
!!    "q_1" : initial right vector
!!
!!    "normalization" : biorthonormalization procedure
!!                      Must be either 'symmetric' or 'asymmetric'
!!
      use array_utilities, only: zero_array
!
      implicit none
!
      type(asymmetric_lanczos_tool) :: lanczos
!
      integer, intent(in) :: n_parameters, chain_length
!
      character (len=200), intent(in) :: normalization
!
      real(dp), dimension(n_parameters), intent(in) :: p1, q1
!
      real(dp) :: ddot, norm
!
      lanczos%n_parameters = n_parameters
      lanczos%chain_length = chain_length
      lanczos%normalization = normalization
!
      call mem%alloc(lanczos%alpha_, lanczos%chain_length)
      call zero_array(lanczos%alpha_, lanczos%chain_length)
!
      call mem%alloc(lanczos%beta_, lanczos%chain_length-1)
      call zero_array(lanczos%beta_, lanczos%chain_length-1)
!      
      call mem%alloc(lanczos%gamma_, lanczos%chain_length-1)
      call zero_array(lanczos%gamma_, lanczos%chain_length-1)
!
      norm=ddot(lanczos%n_parameters, p1, 1, q1, 1)
!
      if (abs(norm - one) .gt. 1.0d-8) &
          call output%error_msg('Start vectors for Asymmetric Lanczos not binormalized')
!      
      lanczos%p_file = sequential_file('asymmetric_lanczos_p', 'unformatted' )
      lanczos%q_file = sequential_file('asymmetric_lanczos_q', 'unformatted' )
      lanczos%T_file = sequential_file('asymmetric_lanczos_T', 'unformatted' )
!
      call lanczos%p_file%open_('write', 'rewind')
      call lanczos%q_file%open_('write', 'rewind')
!
      call lanczos%p_file%write_(p1, lanczos%n_parameters)
      call lanczos%q_file%write_(q1, lanczos%n_parameters)
!
      call lanczos%p_file%close_()
      call lanczos%q_file%close_()
!
   end function
!
!
   subroutine read_p_asymmetric_lanczos_tool(lanczos, p, n)
!!
!!    Read p vector
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!    Adapted for the asymmetric Lanczos solver by Torsha Moitra, Sep 2019
!!
!!    Reads the nth trial vector from file and places it in c.
!!
!!    If n is not passed, it reads the trial at wherever the cursor
!!    in the file currently is.
!!
      implicit none
!
      class(asymmetric_lanczos_tool), intent(inout) :: lanczos
!
      real(dp), dimension(lanczos%n_parameters), intent(out) :: p
!
      integer, optional, intent(in) :: n
!
      call lanczos%p_file%open_('read', 'append')
!
      if (present(n)) then
!
         call lanczos%p_file%rewind_()
         call lanczos%p_file%skip(n - 1)
!
      endif
!
      call lanczos%p_file%read_(p, lanczos%n_parameters)
      call lanczos%p_file%close_()
!
   end subroutine read_p_asymmetric_lanczos_tool


   subroutine read_q_asymmetric_lanczos_tool(lanczos, q, n)
!!
!!    Read q vector
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!    Adapted for the asymmetric Lanczos solver by Torsha Moitra, Sep 2019
!!
!!    Reads the nth trial vector from file and places it in c.
!!
!!    If n is not passed, it reads the trial at wherever the cursor
!!    in the file currently is.
!!
      implicit none
!
      class(asymmetric_lanczos_tool), intent(inout) :: lanczos
!
      real(dp), dimension(lanczos%n_parameters), intent(out) :: q
!
      integer, optional, intent(in) :: n
!
      call lanczos%q_file%open_('read','append')
!
      if (present(n)) then
!
         call lanczos%q_file%rewind_()
         call lanczos%q_file%skip(n - 1)
!
      endif
!
      call lanczos%q_file%read_(q, lanczos%n_parameters)
      call lanczos%q_file%close_()
!
   end subroutine read_q_asymmetric_lanczos_tool

!
   subroutine calculate_alpha_asymmetric_lanczos_tool(lanczos, Aq, i)
!!
!!    Calculate alpha 
!!    Written by Torsha Moitra and Sarai D. Folkestad, Sep 2019
!!
!!    Calculates the diagonal elements of the
!!    tridiagonal matrix T 
!!
!!       alpha_i = p_i^T A q_i = p_i^T (Aq_i)
!!
!!    "Aq" : The linear transform of q_i by some matrix A
!!
     implicit none
!
     class(asymmetric_lanczos_tool), intent(inout) :: lanczos
!
     real(dp), dimension(lanczos%n_parameters), intent(in) :: Aq
     real(dp) :: ddot
!
     integer, intent(in) :: i
!
     real(dp), dimension(:), allocatable :: p
!
     call mem%alloc(p, lanczos%n_parameters)
!
     call lanczos%read_p(p, i)
!
     lanczos%alpha_(i) = ddot(lanczos%n_parameters, p, 1, Aq, 1)
!
     call mem%dealloc(p, lanczos%n_parameters)
!    
   end subroutine calculate_alpha_asymmetric_lanczos_tool
!
!
   subroutine construct_r_asymmetric_lanczos_tool(lanczos, i, Aq, r)
!!
!!    Construct r
!!    Written by Torsha Moitra and Sarai D. Folkestad, Sep 2019
!!
!!    Constructs the r-vector
!!
!!       r_i = A q_i - gamma_i-1 q_i-1 - alpha_i q_i
!!
!!    which is used to calculate gamma and beta, and to 
!!    construct q_i+1
!!
!!    "Aq" : The linear transform of q_i by some matrix A
!!
      use array_utilities, only: copy_and_scale, zero_array
!
      implicit none
!
      class(asymmetric_lanczos_tool), intent(inout) :: lanczos
!
      real(dp), dimension(lanczos%n_parameters), intent(in) :: Aq
      real(dp), dimension(lanczos%n_parameters), intent(out) :: r
!
      integer, intent(in) :: i
!
      real(dp), dimension(:), allocatable :: q
!
!     r_i = Aq_i
!
      call copy_and_scale(one, Aq, r, lanczos%n_parameters)
!
      call mem%alloc(q, lanczos%n_parameters)
!
      call lanczos%read_q(q, i)
!
!     r_i -= alpha_i*q_i
!
      call daxpy(lanczos%n_parameters, -lanczos%alpha_(i), q, 1, r, 1)
!
      if (i .gt. 1) then         
!
!        r -= gamma_i-1 q_i-1
!
         call lanczos%read_q(q, i - 1)
!
         call daxpy(lanczos%n_parameters, -lanczos%gamma_(i - 1), q, 1, r, 1)
!
      endif
!
      call mem%dealloc(q, lanczos%n_parameters) 
!
   end subroutine construct_r_asymmetric_lanczos_tool
!
!
   subroutine construct_s_asymmetric_lanczos_tool(lanczos, i, pA, s)
!!
!!    Construct s
!!    Written by Torsha Moitra and Sarai D. Folkestad, Sep 2019
!!
!!    Constructs the s-vector
!!
!!       s_i = p_i^T A - beta_i-1 p_i-1 - alpha_i p_i
!!
!!    which is used to calculate beta and gamma
!!    and to construct p_i+1     
!!
!!    "pA" : The linear transform of p_i by some matrix A^T
!
      use array_utilities, only: copy_and_scale, zero_array
!
      implicit none
!
      class(asymmetric_lanczos_tool), intent(inout) :: lanczos
!
      real(dp), dimension(lanczos%n_parameters), intent(in) :: pA
      real(dp), dimension(lanczos%n_parameters), intent(out) :: s
!
      integer, intent(in) :: i
!
      real(dp), dimension(:), allocatable :: p
!
!     s_i = p_i^T A 
!
      call copy_and_scale(one, pA, s, lanczos%n_parameters)
!
      call mem%alloc(p, lanczos%n_parameters)
!
      call lanczos%read_p(p, i)
!
!     s_i -= alpha_i p_i 
!
      call daxpy(lanczos%n_parameters, -lanczos%alpha_(i), p, 1, s, 1)      
!
      if (i .gt. 1) then  
!
!        s_i -= beta_i-1 p_i-1
! 
         call lanczos%read_p(p, i - 1)
!
         call daxpy(lanczos%n_parameters, -lanczos%beta_(i - 1), p, 1, s, 1)
!
      endif
!
      call mem%dealloc(p,lanczos%n_parameters)
!
   end subroutine construct_s_asymmetric_lanczos_tool
!
!
   subroutine calculate_beta_gamma_p_q_asymmetric_lanczos_tool(lanczos, i, Aq, pA)
!!
!!    Calculate beta and gamma
!!    Written by Torsha Moitra and Sarai D. Folkestad, Sep 2019
!!
!!    Revised by Torsha Moitra and Sonia Coriani, Nov 2019
!!    Calculates beta and gamma
!!
!!       beta_i  = sqrt(abs(s_i^T r_i))
!!       gamma_i = sign(s_i^T r_i) beta_i,
!!
!!    constructs the new p and q vectors
!!
!!       p_i+1 = s_i/gamma_i
!!       q_i+1 = r_i/beta_i,
!!
!!    and write them to the end of the p and q files.
!!
!!    "pA" : The linear transform of p_i by some matrix A^T
!!
!!    "Aq" : The linear transform of q_i by some matrix A
!!
      implicit none
!
      class(asymmetric_lanczos_tool), intent(inout) :: lanczos
!
      real(dp), dimension(lanczos%n_parameters), intent(in) :: pA, Aq
!
      integer, intent(in) :: i
!
      real(dp) :: ddot
      real(dp), parameter :: s_dot_r_threshold = 1.0d-12
!
      real(dp), dimension(:), allocatable :: r, s
      real(dp) :: s_dot_r
!
      call mem%alloc(r,lanczos%n_parameters)
      call mem%alloc(s,lanczos%n_parameters)
!
      call lanczos%construct_r(i, Aq, r)
      call lanczos%construct_s(i, pA, s)
!
      s_dot_r = ddot(lanczos%n_parameters, s, 1, r, 1)
      lanczos%beta_(i) = dsqrt(dabs(s_dot_r))
      lanczos%gamma_(i)= sign(lanczos%beta_(i), s_dot_r)
!
!     If the overlap of r and s is below the threshold,
!     the chain is terminated at the current iteration
!
      if (dabs(s_dot_r) .lt. s_dot_r_threshold) then
!
          lanczos%chain_length = i
!
      endif
!
!      prepare new q and p vectors 

      call dscal(lanczos%n_parameters, one/lanczos%beta_(i), r, 1)
      call dscal(lanczos%n_parameters, one/lanczos%gamma_(i), s, 1)
!
!     biorthogonalize again previous vectors and binormalize:
!      
      call lanczos%biorthonormalize(r, s, i)
!
      call lanczos%q_file%open_('write', 'append')
      call lanczos%p_file%open_('write', 'append')

      call lanczos%q_file%write_(r, lanczos%n_parameters)
      call lanczos%p_file%write_(s, lanczos%n_parameters)
!    
      call lanczos%p_file%close_()
      call lanczos%q_file%close_()
!
      call mem%dealloc(r,lanczos%n_parameters)
      call mem%dealloc(s,lanczos%n_parameters)
!
   end subroutine calculate_beta_gamma_p_q_asymmetric_lanczos_tool
!
!
   subroutine diagonalize_T_asymmetric_lanczos_tool(lanczos, left, right, &
                  eigenvalues_Re, eigenvalues_Im)
!!
!!    Diagonalize T
!!    Written by Torsha Moitra and Sarai D. Folkestad, Sep 2019
!!
!!    Constructs the tridiagonal matrix T 
!!    which has alpha on the diagonal, beta below the diagonal
!!    and gamma above the diagonal and diagonalize it.
!!
!!    For further information see S. Coriani et al. JCTC, 8, 1616-1628 (2012) 
!!
!!    "left" : On exit contains the left eigenvectors of T
!!
!!    "right" : On exit contains the right eigenvectors of T
!!
!!    "eigenvalues_Re" : On exit contains the real part
!!                       of the eigenvalues of T
!!
!!    "eigenvalues_Im" : On exit contains the imaginary part
!!                       of the eigenvalues of T
!!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(asymmetric_lanczos_tool), intent(inout) :: lanczos
!
      real(dp), dimension(lanczos%chain_length,lanczos%chain_length), intent(out) :: left, right
      real(dp), dimension(lanczos%chain_length), intent(out) :: eigenvalues_Re, eigenvalues_Im
!
      real(dp), dimension(:,:), allocatable :: T
      real(dp), dimension(:), allocatable :: work
!
      integer :: lwork, i, info
!    
      call mem%alloc(T, lanczos%chain_length, lanczos%chain_length)
!
      call zero_array(T, lanczos%chain_length**2)
!
      T(1,1) = lanczos%alpha_(1)
!
!$omp parallel do private(i)
      do i = 2, lanczos%chain_length
!
         T(i, i)     = lanczos%alpha_(i)
         T(i - 1, i) = lanczos%gamma_(i-1)
         T(i, i - 1) = lanczos%beta_(i-1)
!
      enddo
!$omp end parallel do
!
      lwork = 2*lanczos%chain_length**2
!
      call mem%alloc(work, lwork)
!
      call dgeev('V', 'V',                 &
                 lanczos%chain_length, &
                 T,                        &
                 lanczos%chain_length, & 
                 eigenvalues_Re,           &      
                 eigenvalues_Im,           &      
                 left,                     &
                 lanczos%chain_length, &
                 right,                    &
                 lanczos%chain_length, &
                 work,                     &
                 lwork,                    &
                 info)
!
      if (info .ne. 0) then
!
         call output%error_msg('Diagonalization of T ended with error: (i0)', ints=[info])
!
      endif
!
      call mem%dealloc(work,lwork)
      call mem%dealloc(T, lanczos%chain_length, lanczos%chain_length) 
!
   end subroutine diagonalize_t_asymmetric_lanczos_tool
!
!
   subroutine save_T_asymmetric_lanczos_tool(lanczos)
!!
!!    Save alpha, beta, gamma onto T_file
!!    Written by Torsha Moitra and Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(asymmetric_lanczos_tool), intent(inout) :: lanczos
!
      call lanczos%T_file%open_('write', 'rewind')   
!
      call lanczos%T_file%write_(lanczos%alpha_,lanczos%chain_length)
      call lanczos%T_file%write_(lanczos%beta_,lanczos%chain_length-1)
      call lanczos%T_file%write_(lanczos%gamma_,lanczos%chain_length-1)
!
      call lanczos%T_file%close_()

   end subroutine save_T_asymmetric_lanczos_tool
!
!
   subroutine biorthonormalize_asymmetric_lanczos_tool(lanczos,q,p,i)
!!
!!    Biorthonomalize p and q vectors
!!    Written by Torsha Moitra and Sarai D. Folkestad, Sep 2019
!!    Revised by Torsha Moitra and Sonia Coriani, Nov 2019
!!    Enforces the biorthonormality condition on the p and q vectors
!!
!!       P^T*Q = I , 
!!
!!    with 
!!
!!       P = [p_1, p_2, ..., p_i+1], Q = [q_1, q_2, ..., q_i+1].
!!
!!    "q" : Is the newly constructed q_i+1 vector
!!
!!    "p" : Is the newly constructed p_i+1 vector
!
      implicit none
!
      class(asymmetric_lanczos_tool), intent(inout) :: lanczos
!
      real(dp), dimension(lanczos%n_parameters), intent(inout) :: p
      real(dp), dimension(lanczos%n_parameters), intent(inout) :: q
!
      real(dp), dimension(:), allocatable :: p_j, q_j
      real(dp) :: ddot
      real(dp) :: overlap, factor
      
      integer, intent(in) :: i
      integer :: j,k
!      
      call mem%alloc(p_j, lanczos%n_parameters)
      call mem%alloc(q_j, lanczos%n_parameters)
!
      do k = 1,2
!
!       do it twice for numerical stability

        do j = 1, i
!
!         biorthogonalize against previous 

          call lanczos%read_p(p_j, j)      
          call lanczos%read_q(q_j, j)      

          overlap = ddot(lanczos%n_parameters, p_j, 1, q, 1) 
          call daxpy(lanczos%n_parameters, -overlap, q_j, 1, q, 1)   
   
          overlap = ddot(lanczos%n_parameters, p, 1, q_j, 1)
          call daxpy(lanczos%n_parameters, -overlap, p_j, 1, p, 1) 

        enddo

!       now binormalize new p, q pair
!
        if(lanczos%normalization=='symmetric')then
!
          overlap = ddot(lanczos%n_parameters, p, 1, q, 1)
          factor = one/dsqrt(dabs(overlap))
!
          call dscal(lanczos%n_parameters, factor, q, 1)
          call dscal(lanczos%n_parameters, dsign(factor,overlap), p, 1)
!
        else if (lanczos%normalization=='asymmetric')then
!
          overlap = ddot(lanczos%n_parameters, p, 1, q, 1)
          factor= one/dabs(overlap)
          call dscal(lanczos%n_parameters, dsign(factor,overlap), p, 1) 
!
        else 
!
          call output%error_msg('could not recognize biorthonormalization procedure in asymmetric &
            &Lanczos tool.')      
!
        endif
!
      enddo 

      call mem%dealloc(p_j, lanczos%n_parameters)
      call mem%dealloc(q_j, lanczos%n_parameters)
!
   end subroutine biorthonormalize_asymmetric_lanczos_tool
!
!
   subroutine cleanup_asymmetric_lanczos_tool(lanczos)
!!
!!    Cleanup lanczos tool 
!!    Written by Torsha Moitra, Nov 2019
!!
      implicit none
!
      class(asymmetric_lanczos_tool) :: lanczos
      
      call mem%dealloc(lanczos%alpha_, lanczos%chain_length) 
      call mem%dealloc(lanczos%beta_, lanczos%chain_length-1) 
      call mem%dealloc(lanczos%gamma_, lanczos%chain_length-1) 

   end subroutine cleanup_asymmetric_lanczos_tool
!
!
end module asymmetric_lanczos_tool_class
