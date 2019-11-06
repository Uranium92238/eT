!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
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
   module subroutine jacobian_transpose_doubles_a1_doubles_complex(wf, sigma_ai, c_bj, u)
!!
!!    Jacobian transpose doubles A1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad 
!!    and Alexander C. Paul, Feb 2019
!!
      implicit none
!
      class(doubles) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o), intent(in)     :: c_bj
      complex(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u
!
   end subroutine jacobian_transpose_doubles_a1_doubles_complex
!
!
  module subroutine jacobian_transpose_doubles_b1_doubles_complex(wf, sigma_ai, c_bjck)
!!
!!    Jacobian transpose doubles B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad 
!!    and Alexander C. Paul, Feb 2019
!!
      implicit none
!
      class(doubles) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_bjck
      complex(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
!
   end subroutine jacobian_transpose_doubles_b1_doubles_complex
!
!
  module subroutine jacobian_transpose_doubles_a2_doubles_complex(wf, sigma_aibj, c_ai)
!!
!!    Jacobian transpose CC2 A2
!!    Written by Sarai D. Folkestad and Alexander C. Paul, Feb 2019
!!
      implicit none
!
      class(doubles) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: sigma_aibj
!
   end subroutine jacobian_transpose_doubles_a2_doubles_complex
!
