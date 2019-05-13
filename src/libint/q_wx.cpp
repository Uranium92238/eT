//
//
//  eT - a coupled cluster program
//  Copyright (C) 2016-2019 the authors of eT
//
//  eT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  eT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// -----------------------------------------------------------------------
/*
 
  	Construct quadrupole integrals q
  	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2019
 
*/
#include <iostream>

using namespace std;

#include <fstream>
#include <string>
#include <vector>
#include "q_wx.h"
#include "extract_integrals.h"

#include <libint2.hpp>

#include "globals.h"

using namespace libint2;

void construct_ao_q_wx(double *q_xx, double *q_xy, double *q_xz, 
               double *q_yy, double *q_yz, double *q_zz, int *s1, int *s2){

   int thread = omp_get_thread_num();

   const auto& buf_vec = quadrupole[thread].results(); // will point to computed shell sets

   auto n1 = basis[*s1 - 1].size();        // Number of basis functions in shell 1
   auto n2 = basis[*s2 - 1].size();        // number of basis functions in shell 2

   quadrupole[thread].compute(basis[*s1 - 1], basis[*s2 - 1]);

   // Libint computes the overlap integrals and dipole integrals simultaneously
   // We throw these away in this routine
   //
   // buf_vec[0] -> overlap integrals 
   // buf_vec[1] -> mu_x
   // buf_vec[2] -> mu_y
   // buf_vec[3] -> mu_z
   //

   auto ints_shellset_xx = buf_vec[4]; 
   auto ints_shellset_xy = buf_vec[5]; 
   auto ints_shellset_xz = buf_vec[6]; 
   auto ints_shellset_yy = buf_vec[7]; 
   auto ints_shellset_yz = buf_vec[8]; 
   auto ints_shellset_zz = buf_vec[9]; 

   extract_integrals(q_xx, ints_shellset_xx, n1, n2, -1.0e0);
   extract_integrals(q_xy, ints_shellset_xy, n1, n2, -1.0e0);
   extract_integrals(q_xz, ints_shellset_xz, n1, n2, -1.0e0);
   extract_integrals(q_yy, ints_shellset_yy, n1, n2, -1.0e0);
   extract_integrals(q_yz, ints_shellset_yz, n1, n2, -1.0e0);
   extract_integrals(q_zz, ints_shellset_zz, n1, n2, -1.0e0);

   return;

}

