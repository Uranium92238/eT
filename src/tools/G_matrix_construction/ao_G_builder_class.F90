!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module ao_G_builder_class
!
!!
!!    AO G builder class
!!    Written by Sarai D. Folkestad
!!
!!    Directs the building of the two-electron part of the Fock matrix (G)
!!    in the AO basis.
!!
!!    Uses a screening and and adding tool
!!
!
   use parameters
!
   use ao_tool_class,               only: ao_tool
   use ao_eri_getter_class,         only: ao_eri_getter
!
   use memory_manager_class,        only: mem
   use array_initialization,        only: zero_array
   use reordering,                  only: symmetric_sum
!
   use abstract_G_screener_class,   only: abstract_G_screener
   use abstract_G_adder_class,      only: abstract_G_adder
!
   use omp_lib
!
   type :: ao_G_builder
!
      class(abstract_G_screener), allocatable :: screener
      class(abstract_G_adder),   allocatable  :: adder
!
   contains
!
      procedure, private :: construct_threaded
      procedure, public :: construct => construct_ao_G_builder
!
   end type ao_G_builder
!
   interface  ao_G_builder
!
      procedure :: new_ao_G_builder
!
   end interface  ao_G_builder
!
contains
!
!
   pure function new_ao_G_builder(screener, adder) result(this)
!!
!!    New AO G builder
!!    Written by Sarai D. Folkestad, Jun 2021
!!
      implicit none
!
      class(abstract_G_screener), intent(in) :: screener
      class(abstract_G_adder), intent(in) :: adder
      type(ao_G_builder) :: this
!
      this%screener = screener
      this%adder = adder
!
   end function new_ao_G_builder
!
!
   subroutine construct_ao_G_builder(this, ao, eri_getter, D, G)
!!
!!    Construct
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2021
!!
!!    Calculates the two-electron part of the Fock matrix
!!
!!       G(D)_wx = sum_yz g_wxyz D_yz - 1/2 * sum_yz g_wzyx D_yz,
!!
      implicit none
!
      class(ao_G_builder), intent(inout):: this
      class(ao_tool), intent(in):: ao
      class(ao_eri_getter), intent(in) :: eri_getter
!
      real(dp), dimension(ao%n, ao%n), intent(in)    :: D
      real(dp), dimension(ao%n, ao%n), intent(inout) :: G
!
      call this%screener%prescreening(ao, D)
      call this%construct_threaded(eri_getter, ao, D, G)
      call this%screener%cleanup()
!
   end subroutine construct_ao_G_builder
!
!
   subroutine construct_threaded(this, eri_getter, ao, D, G)
!!
!!    Construct threaded
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, and Linda Goletto, 2018-2021
!!
!!    Routine partly based on Hartree-Fock implementation shipped with
!!    the Libint 2 integral package by E. Valeev.
!!
!!    Defines the loop structure to construct the two-electron terms of the Fock matrix
!!
!!        G_wx =+ sum_yz 2g_wxyz D_yz - sum_yz g_wzyx D_yz (= G(D)_wx),
!!
!!    Depending on which version of the adder/screenin tools are used,
!!    either Coulomb, exchange or both are constructed.
!!
!!    Note: contributions from each thread are added to a separate part of G.
!!          Memory scaling is n_threads * n_ao^2
!!
      implicit none
!
      class(ao_G_builder),             intent(in)  :: this
      class(ao_tool),                  intent(in)  :: ao
      class(ao_eri_getter),            intent(in)  :: eri_getter
      real(dp), dimension(ao%n, ao%n), intent(in)  :: D
      real(dp), dimension(ao%n, ao%n), intent(out) :: G
!
      real(dp) :: degeneracy, degeneracy_12, degeneracy_34, degeneracy_12_34
!
      integer :: s1, s2, s3, s4, s4_max, s4_tilde, s1s2, s3s4, s1s2_packed, s3_tilde, s1s2_tilde
!
      real(dp), dimension(ao%max_sh_size**4) :: eri
!
      integer :: skip
!
      real(dp) :: precision_threshold, D_max_s1s2s3s4

      integer :: thread, n_threads
!
      real(dp), dimension(:,:,:), allocatable :: G_thread
!
      precision_threshold = ao%get_libint_epsilon()
!
      n_threads = 1 ! needed when OMP is not enabled
!$    n_threads = omp_get_max_threads()
!
      call mem%alloc(G_thread, ao%n, ao%n, n_threads, set_zero=.true.)
!
!$omp parallel do private(s1, s2, s3, s3_tilde, s4, s1s2, s3s4, s1s2_packed, s4_max, s4_tilde, eri, &
!$omp degeneracy_12, degeneracy_34, degeneracy_12_34, degeneracy, skip, thread, s1s2_tilde, D_max_s1s2s3s4) &
!$omp schedule(dynamic)
      do s1s2_tilde = 1, this%screener%n_sig_s1s2
!
         s1s2 = this%screener%get_sig_s1s2(s1s2_tilde)
!
         s1s2_packed = ao%cs_eri_max_indices(s1s2, 3)
         s1 = ao%cs_eri_max_indices(s1s2_packed, 1)
         s2 = ao%cs_eri_max_indices(s1s2_packed, 2)
!
         thread = 0 ! needed when OMP is not enabled
!$       thread = omp_get_thread_num()
!
         degeneracy_12 = real(2-s2/s1, kind=dp)
!
         do s3_tilde = 1, this%screener%get_n_sig_s3(s1)
!
            s3 = this%screener%get_sig_s3(s3_tilde, s1)
!
            if (this%screener%s1s2s3_exit(ao, s1, s2, s3, s1s2)) exit
            if (this%screener%s1s2s3_cycle(ao, s1, s2, s3, s1s2)) cycle
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4_tilde = 1, ao%n_sig_s2_for_s1(s3)
!
               s4 = ao%sig_s2_for_s1(s4_tilde, s3)
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4
!
               if (s4 .gt. s4_max) exit
               if (this%screener%s1s2s3s4_cycle(ao, s1, s2, s3, s4, s1s2, s3s4)) cycle
!
               D_max_s1s2s3s4 = this%screener%get_D_max(s1, s2, s3, s4)
!
               call eri_getter%get_eri(eri,                            &
                              s1, s2, s3, s4,                          &
                              precision_threshold/(D_max_s1s2s3s4**2), &
                              skip)
!
               if (skip == 1) cycle
!
               degeneracy_34    = real(2-s4/s3, kind=dp)
               degeneracy_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)
               degeneracy       = degeneracy_12*degeneracy_34*degeneracy_12_34
!
               call this%adder%add(ao,                &
                                   D,                 &
                                   G_thread,          &
                                   eri,               &
                                   s1, s2, s3, s4,    &
                                   thread, n_threads, &
                                   degeneracy)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call zero_array(G, (ao%n**2))
      do thread = 1, n_threads
         call daxpy(ao%n**2, one, G_thread(1, 1, thread), 1, G, 1)
      enddo
!
      call mem%dealloc(G_thread, ao%n, ao%n, n_threads)
!
      call symmetric_sum(G, ao%n)
      call dscal(ao%n**2, half, G, 1)
!
   end subroutine construct_threaded
!
!
end module ao_G_builder_class
