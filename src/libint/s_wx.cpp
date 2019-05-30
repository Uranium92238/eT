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
 
  	Overlap matrix routines
  	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
 
*/
#include <iostream>

using namespace std;

#include <fstream>
#include <string>
#include <vector>
#include "s_wx.h"
#include "extract_integrals.h"

#include <libint2.hpp>

#include "globals.h"

using namespace libint2;

void construct_ao_s_wx(double *s, int *s1, int *s2){

   int thread = omp_get_thread_num();

   const auto& buf_vec = overlap[thread].results(); // Will point to computed shell sets

   auto n1 = basis[*s1 - 1].size();          // Number of basis functions in shell 1
   auto n2 = basis[*s2 - 1].size();          // Number of basis functions in shell 2

   overlap[thread].compute(basis[*s1 - 1], basis[*s2 - 1]);

   auto ints_shellset = buf_vec[0]; // location of the computed integrals

   extract_integrals(s, ints_shellset, n1, n2, 1.0e0);

   return;

}

void construct_ao_s_wx_1der(double *s_1x, double *s_1y, double *s_1z, 
                  double *s_2x, double *s_2y, double *s_2z, int *s1, int *s2){

  //
  // This routine constructs the overlap derivatives s_wx^(1) for w in shell s1 and x in shell s2.
  // Here, we have s_wx^(1) = (s_1x, s_1y, ..., s_2z). The elements in s_1x are the derivatives of s_wx
  // with respect to the x-coordinate of the atom that s1 is centered on, etc.
  //

  const auto& buf_vec = overlap_1der.results();

  auto n1 = basis[*s1 - 1].size();
  auto n2 = basis[*s2 - 1].size();

  overlap_1der.compute(basis[*s1 - 1], basis[*s2 - 1]);

  // Get pointers to location of integrals

  auto ints_shellset_1x = buf_vec[0];
  auto ints_shellset_1y = buf_vec[1];
  auto ints_shellset_1z = buf_vec[2];
  auto ints_shellset_2x = buf_vec[3];
  auto ints_shellset_2y = buf_vec[4];
  auto ints_shellset_2z = buf_vec[5];

  // Extract the integrals from each set

  extract_integrals(s_1x, ints_shellset_1x, n1, n2, 1.0e0);
  extract_integrals(s_1y, ints_shellset_1y, n1, n2, 1.0e0);
  extract_integrals(s_1z, ints_shellset_1z, n1, n2, 1.0e0);
  extract_integrals(s_2x, ints_shellset_2x, n1, n2, 1.0e0);
  extract_integrals(s_2y, ints_shellset_2y, n1, n2, 1.0e0);
  extract_integrals(s_2z, ints_shellset_2z, n1, n2, 1.0e0);

  return;

}

