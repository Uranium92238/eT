module kinds
!
!!
!!		Kinds module
!!		Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!  	Defines a set of kinds in terms of precision. Usage:
!!
!! 		In declarations:   "real(dp) :: foo", "integer(i15) :: foo_int", etc.
!! 		In record-lengths: "recl=dp*200" (200 double precision numbers per record)
!!
!
	implicit none
!
	integer, parameter :: sp  = selected_real_kind(6,37)
	integer, parameter :: dp  = selected_real_kind(15,307)
	integer, parameter :: qp  = selected_real_kind(33,4931)
	integer, parameter :: i15 = selected_int_kind(15)
!
end module kinds
