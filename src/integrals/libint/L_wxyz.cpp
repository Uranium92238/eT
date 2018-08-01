/*
/
/ 	Two-electron L_αβγδ routines (for HF)
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
*/
#include <iostream>

using namespace std;

#include <fstream>
#include <string>
#include <vector>
#include <ctime>
#include <libint2.hpp>
#include "globals.h"

#include "utils.h"
#include "L_wxyz.h"

using namespace libint2;

void get_ao_L_wxyz(double *L, long *s1, long *s3){

	const long sh1 = *s1 - 1;
	const long sh3 = *s3 - 1; // C++ arrays start at index 0!

	//initialize();

	int thread = omp_get_thread_num();

	int num_aos = 0;
	num_aos = basis.nbf();

	auto shell2bf = basis.shell2bf(); // maps shell index to basis function index
                                   // shell2bf[0] = index of the first basis function in shell 0
                                   // shell2bf[1] = index of the first basis function in shell 1
                                   // ...

	const auto& buf_vec = electronic_repulsion_engines[thread].results(); // will point to computed shell sets

	const auto bf1 = shell2bf[sh1];  // First basis function in shell 1
	const auto n1 = basis[sh1].size(); // Number of basis functions in shell 1

	const auto bf3 = shell2bf[sh3];  // First basis function in shell 3
	const auto n3 = basis[sh3].size(); // Number of basis functions in shell 3

  	for(auto s2=0; s2!=basis.nshells(); ++s2) {

		auto bf2 = shell2bf[s2];  // First basis function in shell 2
		auto n2 = basis[s2].size(); // Number of basis functions in shell 2

		for (auto s4=0; s4!=basis.nshells(); ++s4) {

			auto bf4 = shell2bf[s4];  // First basis function in shell 4
			auto n4 = basis[s4].size(); // Number of basis functions in shell 4

			electronic_repulsion_engines[thread].compute(basis[sh1], basis[s2], basis[sh3], basis[s4]);

			auto ints_1234 = buf_vec[0]; // Location of computed integrals
    		if (ints_1234 == nullptr)
      					continue;  // nullptr returned if all integrals were screened out

    		for(auto f1=0, f1234=0; f1!=n1; ++f1){
      		for(auto f2=0; f2!=n2; ++f2){

					auto bf2_ind = bf2 + f2;

      			for(auto f3=0; f3!=n3; ++f3){
      				for(auto f4=0; f4!=n4; ++f4, ++f1234){

							auto bf4_ind = bf4 + f4;

				    		*(L - 1 + index_four(f1+1,bf2_ind+1,f3+1,bf4_ind+1,n1,num_aos,n3)) = *(L - 1 + index_four(f1+1,bf2_ind+1,f3+1,bf4_ind+1,n1,num_aos,n3)) + 2.0e0*ints_1234[f1234];
							*(L - 1 + index_four(f1+1,bf4_ind+1,f3+1,bf2_ind+1,n1,num_aos,n3)) = *(L - 1 + index_four(f1+1,bf4_ind+1,f3+1,bf2_ind+1,n1,num_aos,n3)) - 1.0e0*ints_1234[f1234];

						}
					}
				}
			}
		}
	}

	//finalize();

	return;
}
