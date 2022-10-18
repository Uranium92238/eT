!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
module cc_visualization_task_class
!
!!
!! CC visualization task class
!! Written by Eirik F. Kjønstad, 2021
!!
!
   use parameters
   use ccs_class,             only: ccs
   use cc_task_class,         only: cc_task
   use visualization_class,   only: visualization
   use memory_manager_class,  only: mem
!
   implicit none
!
   type, extends(cc_task) :: cc_visualization_task
!
      logical, private  :: plot_gs_density, plot_mn_densities, plot_es_densities
      logical, private  :: plot_ntos, plot_cntos, plot_mos
!
      real(dp), private :: nto_threshold
!
      character(len=5) :: es_transformation
!
      integer :: n_states_to_plot_densities, n_states_to_plot_ntos
      integer, dimension(:), allocatable :: states_to_plot_densities, states_to_plot_ntos
!
      class(visualization), allocatable, private :: visualizer
!
   contains
!
      procedure, public :: execute &
                        => execute_cc_visualization_task
!
      procedure, private :: set_states_to_plot_densities
      procedure, private :: set_states_to_plot_ntos
      procedure, private :: visualize_ntos
      procedure, private :: plot_orbitals
      procedure, private :: visualize_molecular_orbitals
      procedure, private :: visualize_ground_state_density
      procedure, private :: visualize_transition_densities
      procedure, private :: visualize_excited_state_densities
!
      procedure, private :: cleanup
!
   end type cc_visualization_task
!
!
   interface cc_visualization_task
!
      procedure :: new_cc_visualization_task
!
   end interface cc_visualization_task
!
!
contains
!
!
   function new_cc_visualization_task() result(this)
!!
!!    New CC visualization task
!!    Written by Eirik F. Kjønstad, 2021
!!
      use global_in, only: input
!
      implicit none
!
      type(cc_visualization_task) :: this
!
      this%name_ = 'Plotting orbitals and/or CC densities'
!
      this%plot_mos = input%is_keyword_present('plot cc orbitals', 'visualization')
!
      this%plot_gs_density = input%is_keyword_present('plot cc density', 'visualization')
!
      this%plot_mn_densities = input%is_keyword_present('plot transition densities', &
                                                        'visualization')
!
      this%plot_es_densities = input%is_keyword_present('plot es densities', &
                                                        'visualization')
!
      this%plot_ntos  = input%is_keyword_present('plot ntos','visualization')
!
      this%plot_cntos = input%is_keyword_present('plot cntos','visualization')
!
      this%nto_threshold = 0.1d0
      call input%get_keyword('nto threshold', 'visualization', this%nto_threshold)
!
      this%es_transformation = "right"
!
      if (input%is_keyword_present('left eigenvectors', 'solver cc es')) then
         this%es_transformation = "left"
      end if
!
      if (input%is_keyword_present('right eigenvectors', 'solver cc es')) then
         this%es_transformation = "right"
      end if
!
   end function new_cc_visualization_task
!
!
   subroutine execute_cc_visualization_task(this, wf)
!!
!!    Execute
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(cc_visualization_task), intent(inout) :: this
!
      class(ccs), intent(inout), target :: wf
!
      logical :: plot_anything
      character(len=4) :: tag
!
      plot_anything = any([this%plot_mos, this%plot_gs_density, this%plot_mn_densities, &
                           this%plot_es_densities, this%plot_ntos, this%plot_cntos])
!
!
      if (.not. plot_anything) return
!
      call this%print_header()
      call this%start_timer()
!
      call this%set_states_to_plot_densities()
!
      this%visualizer = visualization(wf%ao)
      call this%visualizer%initialize(wf%ao)
!
      call this%visualize_molecular_orbitals(wf)
!
      call this%visualize_ground_state_density(wf)
!
      call this%visualize_transition_densities(wf)
      call this%visualize_excited_state_densities(wf)
!
      if (this%plot_ntos .or. this%plot_cntos) then
!
         if (this%plot_ntos) tag = 'nto'
         if (this%plot_cntos) tag = 'cnto'
         call this%set_states_to_plot_ntos(trim(tag), wf%n_singlet_states)
         call this%visualize_ntos(wf, trim(tag))
!
      end if
!
      call this%visualizer%cleanup()
!
      call this%cleanup()
!
      call this%end_timer()
!
   end subroutine execute_cc_visualization_task
!
!
   subroutine visualize_molecular_orbitals(this, wf)
!!
!!    Visualize molecular orbitals
!!    Written by Sarai D. Folkestad and Alexander C. Paul, 2019-2022
!!
      use global_in, only: input
!
      implicit none
!
      class(cc_visualization_task), intent(inout) :: this
      class(ccs), intent(in) :: wf
!
      integer :: i, n_orbitals_to_plot
!
      integer, dimension(:), allocatable :: orbitals_to_plot
!
      real(dp), dimension(:,:), allocatable :: coefficients
!
      character(len=200), dimension(:), allocatable :: orbital_file_tags
!
      if (.not. this%plot_mos) return
!
      n_orbitals_to_plot = input%get_n_elements_for_keyword('plot cc orbitals', 'visualization')
!
      call mem%alloc(orbitals_to_plot, n_orbitals_to_plot)
      call input%get_array_for_keyword('plot cc orbitals', 'visualization', &
                                        n_orbitals_to_plot, orbitals_to_plot)
!
      call mem%alloc(coefficients, wf%ao%n, n_orbitals_to_plot)
      allocate(orbital_file_tags(n_orbitals_to_plot))
!
      do i = 1, n_orbitals_to_plot
!
         call dcopy(wf%ao%n, wf%orbital_coefficients(1, orbitals_to_plot(i)), 1, &
                    coefficients(1, i), 1)
!
         write(orbital_file_tags(i), '(a, i4.4)') 'CC_MO_', orbitals_to_plot(i)
!
      enddo
!
      call mem%dealloc(orbitals_to_plot, n_orbitals_to_plot)
!
      call this%visualizer%plot_orbitals(wf%ao, coefficients, &
                                         n_orbitals_to_plot, orbital_file_tags)
!
      call mem%dealloc(coefficients, wf%ao%n, n_orbitals_to_plot)
      deallocate(orbital_file_tags)
!
   end subroutine visualize_molecular_orbitals
!
!
   subroutine visualize_ground_state_density(this, wf)
!!
!!    Visualize ground state density
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(cc_visualization_task), intent(inout) :: this
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:), allocatable :: c_D_ct
!
      if (.not. this%plot_gs_density) return
!
!     Back-transform density matrix to the AO basis
!
      call mem%alloc(c_D_ct, wf%ao%n, wf%ao%n)
      call wf%add_t1_terms_and_transform(wf%density, c_D_ct)
!
      call this%visualizer%plot_density(wf%ao, c_D_ct, 'cc_gs_density')
!
      call mem%dealloc(c_D_ct, wf%ao%n, wf%ao%n)
!
   end subroutine visualize_ground_state_density
!
!
   subroutine visualize_ntos(this, wf, tag)
!!
!!    Visualize ntos
!!    Written by Sarai D. Folkestad, May 2020
!!
      use global_out, only: output
!
      implicit none
!
      class(cc_visualization_task), intent(inout) :: this
      class(ccs), intent(in) :: wf
      character(len=*), intent(in) :: tag
!
      integer :: state, n_significant_v, n_significant_o, p
!
      real(dp), dimension(:,:), allocatable :: orbitals
!
      if ((this%nto_threshold .gt. 1) .or. (this%nto_threshold .lt. 0)) then
            call output%error_msg('illegal threshold given for NTO/CNTO plotting, (e8.4)', &
                                  reals=[this%nto_threshold])
      end if
!
      call mem%alloc(orbitals, wf%ao%n, wf%n_mo)
!
      do p = 1, this%n_states_to_plot_ntos
!
         state = this%states_to_plot_ntos(p)
!
         call wf%construct_ntos_or_cntos(orbitals, state, n_significant_v, n_significant_o, &
                                         this%es_transformation, tag, this%nto_threshold)
!
         call this%plot_orbitals(wf, n_significant_o, n_significant_v, orbitals, tag, state)
!
      enddo
!
      call mem%dealloc(orbitals, wf%ao%n, wf%n_mo)
!
   end subroutine visualize_ntos
!
!
   subroutine plot_orbitals(this, wf, n_sig_o, n_sig_v, orbitals, tag, state)
!!
!!    Plot orbitals
!!    Written by Sarai D. Folkestad, May 2020
!!
      implicit none
!
      class(cc_visualization_task),          intent(inout) :: this
      class(ccs),                            intent(in) :: wf
!
      integer,                               intent(in) :: n_sig_o, n_sig_v, state
      real(dp), dimension(wf%ao%n, wf%n_mo), intent(in) :: orbitals
      character(len=*),                      intent(in) :: tag
!
      character(len=200), dimension(:), allocatable :: file_tags
!
      integer :: i
      character(len=2) :: l_or_r
!
      l_or_r = this%es_transformation(1:1) // '_'
!
      allocate(file_tags(max(n_sig_o, n_sig_v)))
!
      do i = 1, n_sig_o
!
         write(file_tags(i), '(a,i3.3,3a,i3.3)') l_or_r, state, '_', tag, '_o_', i
!
      enddo
!
      call this%visualizer%plot_orbitals(wf%ao, orbitals(:, 1:n_sig_o), &
                                         n_sig_o, file_tags)
!
      do i = 1, n_sig_v
!
         write(file_tags(i), '(a,i3.3,3a,i3.3)') l_or_r, state, '_', tag, '_v_', i
!
      enddo
!
      call this%visualizer%plot_orbitals(wf%ao, orbitals(:, wf%n_o+1:n_sig_v), &
                                         n_sig_v, file_tags)
!
      deallocate(file_tags)
!
   end subroutine plot_orbitals
!
!
   subroutine visualize_transition_densities(this, wf)
!!
!!    Visualize transition densities
!!    Written by Alexander C. Paul, Dec 2020
!!
      implicit none
!
      class(cc_visualization_task), intent(inout) :: this
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:), allocatable :: c_D_ct, density
!
      integer :: p, state_p, state_q
      logical :: file_read
      character(len=200) :: tag
!
      if (.not. this%plot_mn_densities) return
!
      call mem%alloc(density, wf%n_mo, wf%n_mo)
      call mem%alloc(c_D_ct, wf%ao%n, wf%ao%n)
!
      do p = 1, this%n_states_to_plot_densities
!
         state_p = this%states_to_plot_densities(p)
!
         do state_q = 1, wf%n_singlet_states
!
            if (state_p == state_q) cycle
!
            call wf%get_density_for_plotting(c_D_ct, density, state_p, state_q, file_read)
!
            if (file_read) then
               write(tag, '(a, i3.3, a, i3.3)') 'dm_', state_p, '_', state_q
               call this%visualizer%plot_density(wf%ao, c_D_ct, tag)
            end if
!
            call wf%get_density_for_plotting(c_D_ct, density, state_q, state_p, file_read)
!
            if (file_read) then
               write(tag, '(a, i3.3, a, i3.3)') 'dm_', state_q, '_', state_p
               call this%visualizer%plot_density(wf%ao, c_D_ct, tag)
            end if
!
         end do
      end do
!
      call mem%dealloc(density, wf%n_mo, wf%n_mo)
      call mem%dealloc(c_D_ct, wf%ao%n, wf%ao%n)
!
   end subroutine visualize_transition_densities
!
!
   subroutine visualize_excited_state_densities(this, wf)
!!
!!    Visualize excited state densities
!!    Written by Alexander C. Paul, Dec 2020
!!
      implicit none
!
      class(cc_visualization_task), intent(inout) :: this
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:), allocatable :: c_D_ct, density
!
      integer :: p, state_p
      logical :: file_read
      character(len=200) :: tag
!
      if (.not. this%plot_es_densities) return
!
      call mem%alloc(density, wf%n_mo, wf%n_mo)
      call mem%alloc(c_D_ct, wf%ao%n, wf%ao%n)
!
      do p = 2, this%n_states_to_plot_densities ! first initial state is the GS
!
         state_p = this%states_to_plot_densities(p)
!
         call wf%get_density_for_plotting(c_D_ct, density, state_p, state_p, file_read)
!
         if (file_read) then
            write(tag, '(a, i3.3, a, i3.3)') 'dm_', state_p, '_', state_p
            call this%visualizer%plot_density(wf%ao, c_D_ct, tag)
         end if
!
      end do
!
      call mem%dealloc(density, wf%n_mo, wf%n_mo)
      call mem%dealloc(c_D_ct, wf%ao%n, wf%ao%n)
!
   end subroutine visualize_excited_state_densities
!
!
   subroutine set_states_to_plot_densities(this)
!!
!!    Set states to plot densities
!!    Written Alexander C. Paul, Apr 2021
!!
!!    Determine the states for which densities shall be plotted.
!!    If not present the densities for all initial states are requested
!!
      use global_in, only: input
!
      implicit none
!
      class(cc_visualization_task), intent(inout) :: this
!
      if (input%is_keyword_present('states to plot', 'visualization')) then
!
         this%n_states_to_plot_densities = input%get_n_elements_for_keyword('states to plot', &
                                                                            'visualization')
!
         call mem%alloc(this%states_to_plot_densities, this%n_states_to_plot_densities)
!
         call input%get_array_for_keyword('states to plot', 'visualization', &
                                          this%n_states_to_plot_densities, &
                                          this%states_to_plot_densities)

!
         if (any(this%states_to_plot_densities == 0)) this%plot_gs_density = .true.
!
      else if (input%is_keyword_present('initial states','cc response')) then
!
         this%n_states_to_plot_densities = input%get_n_elements_for_keyword('initial states', &
                                                                            'cc response')
!
         call mem%alloc(this%states_to_plot_densities, this%n_states_to_plot_densities)
!
         call input%get_array_for_keyword('initial states', 'cc response', &
                                          this%n_states_to_plot_densities, &
                                          this%states_to_plot_densities)
!
         if (any(this%states_to_plot_densities == 0)) this%plot_gs_density = .true.
!
      end if
!
   end subroutine set_states_to_plot_densities
!
!
   subroutine set_states_to_plot_ntos(this, tag, default_n_states)
!!
!!    Set states to plot
!!    Written Alexander C. Paul, Apr 2021
!!
!!    Sets number of states for which to plot ntos/cntos
!!    Default: plot ntos for all orbitals
!!
!!    tag: either nto or cnto
!!
      use global_in, only: input
!
      implicit none
!
      class(cc_visualization_task), intent(inout) :: this
      character(len=*), intent(in) :: tag
      integer, intent(in) :: default_n_states
!
      integer :: i
      character(len=:), allocatable :: keyword
!
      keyword = 'plot ' // trim(tag) // 's'
!
      this%n_states_to_plot_ntos = input%get_n_elements_for_keyword(keyword, 'visualization')
!
      if (this%n_states_to_plot_ntos == 0) then
!
         this%n_states_to_plot_ntos = default_n_states
         call mem%alloc(this%states_to_plot_ntos, this%n_states_to_plot_ntos)
         this%states_to_plot_ntos = [(i, i=1, default_n_states)]
         return
!
      end if
!
      call mem%alloc(this%states_to_plot_ntos, this%n_states_to_plot_ntos)
!
      call input%get_array_for_keyword(keyword, 'visualization', &
                                       this%n_states_to_plot_ntos, &
                                       this%states_to_plot_ntos)
!
   end subroutine set_states_to_plot_ntos
!
!
   subroutine cleanup(this)
!!
!!    Cleanup
!!    Written Alexander C. Paul, Apr 2021
!!
      implicit none
!
      class(cc_visualization_task), intent(inout) :: this
!
      if (allocated(this%states_to_plot_densities)) &
         call mem%dealloc(this%states_to_plot_densities, this%n_states_to_plot_densities)
!
      if (allocated(this%states_to_plot_ntos)) &
         call mem%dealloc(this%states_to_plot_ntos, this%n_states_to_plot_ntos)
!
   end subroutine cleanup
!
!
end module cc_visualization_task_class
