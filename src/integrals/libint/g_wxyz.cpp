/*
/
/ 	Two-electron integral routines (for g_wxyz = g_αβγδ)
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
*/
#include <iostream>

using namespace std;

#include <fstream>
#include <string>
#include <ctime>

#include "utils.h"
#include "g_wxyz.h"

#include <libint2.hpp>
#include <vector>

#include "globals.h"

using namespace libint2;

void get_ao_g_wxyz(double *g, long *s1, long *s2, long *s3, long *s4){

	initialize();

	int thread = omp_get_thread_num();

	const long shell1 = *s1 - 1;
	const long shell2 = *s2 - 1;
	const long shell3 = *s3 - 1;
	const long shell4 = *s4 - 1;

	int num_aos = 0;
	num_aos = basis.nbf();

	auto shell2bf = basis.shell2bf(); // maps shell index to basis function index
                                   // shell2bf[0] = index of the first basis function in shell 0
                                   // shell2bf[1] = index of the first basis function in shell 1
                                   // ...

	const auto& buf_vec = electronic_repulsion_engines[thread].results(); // will point to computed shell sets

	clock_t t;
	double tacc;
	tacc = 0;

	auto bf1 = shell2bf[shell1];  // First basis function in shell 1
	auto n1 = basis[shell1].size(); // Number of basis functions in shell 1

	auto bf2 = shell2bf[shell2];  // First basis function in shell 2
	auto n2 = basis[shell2].size(); // Number of basis functions in shell 2

	auto bf3 = shell2bf[shell3];  // First basis function in shell 3
	auto n3 = basis[shell3].size(); // Number of basis functions in shell 3

	auto bf4 = shell2bf[shell4];  // First basis function in shell 4
	auto n4 = basis[shell4].size(); // Number of basis functions in shell 4

	electronic_repulsion_engines[thread].compute(basis[shell1], basis[shell2], basis[shell3], basis[shell4]);

	auto ints_1234 = buf_vec[0]; // Location of computed integrals

   if (ints_1234 == nullptr)
   {
      for(auto f1=0, f1234=0; f1!=n1; ++f1){

  	      for(auto f2=0; f2!=n2; ++f2){

            for(auto f3=0; f3!=n3; ++f3){

               for(auto f4=0; f4!=n4; ++f4, ++f1234){

	  				  *(g - 1 + index_four(f1+1,f2+1,f3+1,f4+1,n1,n2,n3)) = 0.0e0;

               }
            }
         }
      }
   }
   else
   {
      for(auto f1=0, f1234=0; f1!=n1; ++f1){

         for(auto f2=0; f2!=n2; ++f2){

            for(auto f3=0; f3!=n3; ++f3){

               for(auto f4=0; f4!=n4; ++f4, ++f1234){

                 *(g - 1 + index_four(f1+1,f2+1,f3+1,f4+1,n1,n2,n3)) = ints_1234[f1234];

               }
            }
         }
      }
   }

	finalize();

	return;
}
