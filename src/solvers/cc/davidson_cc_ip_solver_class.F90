module davidson_cc_ip_solver_class
!
!!
!!    Davidson coupled cluster ionized state solver class module
!!    Written by Sarai D. Folkestad, 2018
!!
!!    
!
   use davidson_cc_es_solver_class
!
   implicit none
!
   type, extends(davidson_cc_es_solver) :: davidson_cc_ip_solver
!
   contains
!
    ! procedure :: print_banner           => print_banner_davidson_cvs_cc_es_solver
    ! procedure :: print_summary  => print_summary_davidson_cvs_cc_es_solver
!
    ! procedure :: read_settings          => read_settings_davidson_cvs_cc_es_solver
!
    ! procedure :: print_settings => print_settings_davidson_cvs_cc_es_solver
!
    ! procedure :: set_start_vectors      => set_start_vectors_davidson_cvs_cc_es_solver
    ! procedure :: set_projection_vector  => set_projection_vector_davidson_cvs_cc_es_solver
!
    ! procedure :: initialize_core_MOs    => initialize_core_MOs_davidson_cvs_cc_es_solver
    ! procedure :: initialize_cores       => initialize_cores_davidson_cvs_cc_es_solver
!
    ! procedure :: destruct_core_MOs      => destruct_core_MOs_davidson_cvs_cc_es_solver
    ! procedure :: destruct_cores         => destruct_cores_davidson_cvs_cc_es_solver
!
   end type davidson_cc_ip_solver
!
!
contains
!
!
!
!
end module davidson_cc_ip_solver_class
