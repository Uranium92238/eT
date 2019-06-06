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
  
   	Construct h
   	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
  
*/

#include <iostream>

using namespace std;

#include <fstream>
#include <string>
#include <vector>

#include <libint2.hpp>

#include "h_wx.h"
#include "extract_integrals.h"

#include "globals.h"
#include "omp_control.h"

using namespace libint2;

void construct_ao_h_wx(double *h, int *s1, int *s2){
/*
/   Add kinetic contribution
*/
  int thread = omp_get_thread_num();

  const auto& buf_vec = kinetic[thread].results(); // will point to computed shell sets

  kinetic[thread].compute(basis[*s1 - 1], basis[*s2 - 1]);

  auto ints_shellset = buf_vec[0];                // location of the computed integrals

  auto n1 = basis[*s1 - 1].size();                // Number of basis functions in shell 1
  auto n2 = basis[*s2 - 1].size();                // number of basis functions in shell 2

  extract_integrals(h, ints_shellset, n1, n2, 1.0e0);

/*
/   Add nuclear attraction contribution
*/
  const auto& buf_vec_n = nuclear[thread].results();      // will point to computed shell sets

  nuclear[thread].compute(basis[*s1 - 1], basis[*s2 - 1]);
  ints_shellset = buf_vec_n[0];                          // location of the computed integrals

  extract_and_add_integrals(h, ints_shellset, n1, n2, 1.0e0);

  return;
}

void construct_ao_kinetic_1der(double *h_1x, double *h_1y, double *h_1z, 
                  double *h_2x, double *h_2y, double *h_2z, int *s1, int *s2){
/*
    Add kinetic contribution to derivative of h 
*/

  const auto& buf_vec = kinetic_1der.results();

  auto n1 = basis[*s1 - 1].size();
  auto n2 = basis[*s2 - 1].size();

  kinetic_1der.compute(basis[*s1 - 1], basis[*s2 - 1]);

  // Get pointers to location of integrals

  auto ints_shellset_1x = buf_vec[0];
  auto ints_shellset_1y = buf_vec[1];
  auto ints_shellset_1z = buf_vec[2];
  auto ints_shellset_2x = buf_vec[3];
  auto ints_shellset_2y = buf_vec[4];
  auto ints_shellset_2z = buf_vec[5];

  // Extract the integrals from each set

  extract_integrals(h_1x, ints_shellset_1x, n1, n2, 1.0e0);
  extract_integrals(h_1y, ints_shellset_1y, n1, n2, 1.0e0);
  extract_integrals(h_1z, ints_shellset_1z, n1, n2, 1.0e0);
  extract_integrals(h_2x, ints_shellset_2x, n1, n2, 1.0e0);
  extract_integrals(h_2y, ints_shellset_2y, n1, n2, 1.0e0);
  extract_integrals(h_2z, ints_shellset_2z, n1, n2, 1.0e0);

}

void construct_ao_h_wx_1der(double *h_wxqk, int *s1, int *s2){
/*
/   Add nuclear 1st derivative contribution to h
*/
  const auto& ints = nuclear_1der.results(); // will point to computed shell sets

  nuclear_1der.compute(basis[*s1 - 1], basis[*s2 - 1]);

  auto n_atoms = atoms.size();

  auto atom1 = shell2atom[*s1];
  auto atom2 = shell2atom[*s2];
  auto n_ao = basis.size();

  auto n_centers = n_atoms + 2; 

  for (auto center = 0, auto shellset = 0; center != n_centers; ++center){

    auto atom = (center == 0) ? atom1 : ((center == 1) ? atom2 : center - 2);

    for (auto coordinate = 0; coordinate != 3; ++coordinate, ++shellset){

      auto hwqx_offset = n_ao*(n_ao*(3*(atom + 1) + (coordinate + 1))+1)+1;

      

    }

  }

}
