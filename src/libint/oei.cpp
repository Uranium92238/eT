//
//
//  eT - a coupled cluster program
//  Copyright (C) 2016-2021 the authors of eT
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
 
   One-electron integral (oei) routines 
   Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2020
 
*/
#include <iostream>

using namespace std;

#include "oei.h"
#include "extract_integrals.h"
#include "globals.h"


void get_oei(char *type_, double *x, int *s1, int *s2){

/*
   Get one-electron integrals (oei)
   Written by Sarai D. Folkestad, Sep 2020

   Wrapper for the construction of oeis 
   The string 'type_' is used to determine the libint engine
   passed to the 'construct_oei' routine.

*/
   string engine_type   = string(type_); // Convert to string for ease of comparisons
   bool add_            = false;         // Adding to x or copying into x
   double prefactor     = 1.0e0;         // Prefactor for integrals
   int offset           = 0;             // Offset to get the right operator from the engine
   int n_components     = 1;             // Total number of components from the engine: 
                                         //  components operator 1 + components operator 2 + ...

   if (engine_type == "overlap") {

      construct_oei(overlap,        \
                    x,              \
                    s1,             \
                    s2,             \
                    n_components,   \
                    offset,         \
                    prefactor,      \
                    add_);

   }
   else if (engine_type == "hamiltonian") {

      // Set x to kinetic
      construct_oei(kinetic,        \
                    x,              \
                    s1,             \
                    s2,             \
                    n_components,   \
                    offset,         \
                    prefactor,      \
                    add_);

      // Add nuclear to x
      add_ = true;
      construct_oei(nuclear,        \
                    x,              \
                    s1,             \
                    s2,             \
                    n_components,   \
                    offset,         \
                    prefactor,      \
                    add_);
   }
   else if (engine_type == "dipole") {

      prefactor      = -1.0e0;   // eT sign convention
      n_components   = 4;        // n_components = 1 + 3 (overlap + dipole)
      offset         = 1;        // Skip overlap
      construct_oei(dipole,         \
                    x,              \
                    s1,             \
                    s2,             \
                    n_components,   \
                    offset,         \
                    prefactor,      \
                    add_);
   }
   else if (engine_type == "quadrupole") {


      prefactor      = -1.0e0;   // eT sign convention
      n_components   = 10;       // n_components = 1 + 3 + 6 (overlap + dipole + quadrupole)
      offset         = 4;        // Skip overlap and dipole
      construct_oei(quadrupole,     \
                     x,             \
                     s1,            \
                     s2,            \
                     n_components,  \
                     offset,        \
                     prefactor,     \
                     add_);
   }
   else if (engine_type == "electrostatic potential") {

      construct_oei(coulomb_external_charges,   \
                    x,                          \
                    s1,                         \
                    s2,                         \
                    n_components,               \
                    offset,                     \
                    prefactor,                  \
                    add_);

   }
   else if (engine_type == "electrostatic potential unit") {

      construct_oei(coulomb_external_unit_charges, \
                    x,                             \
                    s1,                            \
                    s2,                            \
                    n_components,                  \
                    offset,                        \
                    prefactor,                     \
                    add_);

   }
   return;
}

void get_oei_1der(char *type_, double *x, int *s1, int *s2){

/*
   Get one-electron integrals (oei) 1 der 
   Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2020

   Wrapper for the construction of oeis 
   The string 'type_' is used to determine the libint engine
   passed to the 'construct_oei' routine.

*/
   string engine_type   = string(type_); // Convert to string for ease of comparisons
   bool add_            = false;         // Adding to x or copying into x
   double prefactor     = 1.0e0;         // Prefactor for integrals
   int offset           = 0;             // Offset to get the right operator from the engine
   int n_components     = 6;             // Total number of components from the engine: 
                                         //  components operator 1 + components operator 2 + ...

   if (engine_type == "overlap") {

      construct_oei(overlap_1der,   \
                    x,              \
                    s1,             \
                    s2,             \
                    n_components,   \
                    offset,         \
                    prefactor,      \
                    add_);

   }
   else if (engine_type == "kinetic") {

      construct_oei(kinetic_1der,   \
                    x,              \
                    s1,             \
                    s2,             \
                    n_components,   \
                    offset,         \
                    prefactor,      \
                    add_);

   }
   return;
}

void construct_oei(vector<libint2::Engine> engine, \
                  double *x,                       \
                  int *s1,                         \
                  int *s2,                         \
                  int n_components,                \
                  int offset,                      \
                  double prefactor,                \
                  bool add_){

/*
   Construct one-electron integral(oei)
   Written by Eirik F. Kjønstad and Sarai D. Folkestad 2018-2020

   Constructs the oeis by using the engine passed as an argument.

   x              : Array to hold calculated integrals
   s1, s2         : Shells to compute integrals for
   n_components   : Number of components
   offset         : Offset used to skip additional operators
   prefactor      : Multiplicative prefactor
   add_           : Add or overwrite x
*/

   int thread = omp_get_thread_num();

   const auto& buf_vec = engine[thread].results(); // Will point to computed shell sets

   auto n1 = basis[*s1 - 1].size();          // Number of basis functions in shell 1
   auto n2 = basis[*s2 - 1].size();          // Number of basis functions in shell 2

   engine[thread].compute(basis[*s1 - 1], basis[*s2 - 1]);

   /*
      Loop over the total number of components for the given engine
      and extract those corresponding to the wanted property (given by offset)
   */

   for (int i = offset; i < n_components; i++){

      auto ints_shellset = buf_vec[i]; // location of the computed integrals
      int x_offset = (i - offset)*(n1*n2);

      if (add_) {

         extract_and_add_integrals(&x[x_offset], ints_shellset, n1, n2, prefactor);
      }
      else{

         extract_integrals(&x[x_offset], ints_shellset, n1, n2, prefactor);
      }
   }
   return;
}

void add_nuclear_h_1der(double *h_wxqk, int *s1, int *s2, int *n_ao){
/*

    Add nuclear h 1der
    Written by Eirik F. Kjønstad, 2019-2020

    Adds the contributions to the nuclear first derivative of h to 
    the array h_wxqk (indices: ao, ao, xyz, atom) associated with the shells
    s1 and s2. 

    Note that since the operator depends on the atomic centers, there are 
    contributions to all atoms (i.e. k = 1,2,3,...,n_atoms) stemming from 
    s1 and s2. This is unlike for center-independent operators, see the 
    get_oei_1der routine.

*/

  int thread = omp_get_thread_num();

  const auto& integrals = nuclear_1der[thread].results(); // will point to computed shell sets

  nuclear_1der[thread].compute(basis[*s1 - 1], basis[*s2 - 1]);

  auto shell_to_atom = basis.shell2atom(atoms);

  int atom1 = shell_to_atom[*s1 - 1];
  int atom2 = shell_to_atom[*s2 - 1];

  auto shell2bf = basis.shell2bf();

  int s1_first = shell2bf[*s1-1];
  int s2_first = shell2bf[*s2-1];

  int n1 = basis[*s1 - 1].size();
  int n2 = basis[*s2 - 1].size();

  const auto n_atoms = atoms.size();
  const auto n_centers = n_atoms + 2;

/*

    Libint order 'integrals' shellsets:

    Shellset 1-3: xyz (orbital) derivative contributions with respect to center given by s1 
    Shellset 4-6: xyz (orbital) derivative contributions with respect to center given by s2 
    Shellset 7-:  xyz (operator) derivative contributions due to atoms 1, 2, 3, ..., n_atoms


*/

  int nao = *n_ao;

  for (std::size_t center = 0, shellset = 0; center != n_centers; ++center){

    int atom = (center == 0) ? atom1 : ((center == 1) ? atom2 : center - 2);

    for (int xyz = 0; xyz != 3; ++xyz, ++shellset){

      if (integrals[shellset] != nullptr){

        const auto current_integrals = integrals[shellset];

        for (int f2 = 0; f2 != n2; ++f2){
          for (int f1 = 0; f1 != n1; ++f1){

            auto offset = nao*(nao*(3*(atom)+xyz)+s2_first+f2)+s1_first+f1;
            h_wxqk[offset] = h_wxqk[offset] + current_integrals[f1*n2 + f2];

          }
        }
        if (*s1 != *s2){
          for (int f1 = 0, f12 = 0; f1 != n1; ++f1){
            for (int f2 = 0; f2 != n2; ++f2, ++f12){

               auto offset = nao*(nao*(3*(atom)+xyz)+s1_first+f1)+s2_first+f2;
               h_wxqk[offset] = h_wxqk[offset] + current_integrals[f12];

            }
          }          
        }
      }
    }
  }
  return;
}
