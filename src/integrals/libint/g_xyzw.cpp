/*
/
/ 	Two-electron integral routines (for g_xyzw = g_αβγδ)
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
*/
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>

#include "utils.h"
#include "g_xyzw.h"

#include <libint2.hpp>

#include "globals.h"

using namespace std;
using namespace libint2;

void get_ao_g_xyzw(double *g){

	initialize();

	int num_aos = 0;
	num_aos = basis.nbf();

	auto shell2bf = basis.shell2bf(); // maps shell index to basis function index
                                   // shell2bf[0] = index of the first basis function in shell 0
                                   // shell2bf[1] = index of the first basis function in shell 1
                                   // ...

	const auto& buf_vec = electronic_repulsion.results(); // will point to computed shell sets

	clock_t t;
	double tacc;
	tacc = 0;

	for(auto s1=0; s1!=basis.size(); ++s1) {

		auto bf1 = shell2bf[s1];  // First basis function in shell 1
		auto n1 = basis[s1].size(); // Number of basis functions in shell 1

  		for(auto s2=0; s2!=basis.size(); ++s2) {

			auto bf2 = shell2bf[s2];  // First basis function in shell 2
			auto n2 = basis[s2].size(); // Number of basis functions in shell 2

			for (auto s3=0; s3!=basis.size(); ++s3) {

				auto bf3 = shell2bf[s3];  // First basis function in shell 3
				auto n3 = basis[s3].size(); // Number of basis functions in shell 3

				for (auto s4=0; s4!=basis.size(); ++s4) {

					auto bf4 = shell2bf[s4];  // First basis function in shell 4
					auto n4 = basis[s4].size(); // Number of basis functions in shell 4

					electronic_repulsion.compute(basis[s1], basis[s2], basis[s3], basis[s4]);

					auto ints_1234 = buf_vec[0]; // Location of computed integrals
    				if (ints_1234 == nullptr)
      					continue;  // nullptr returned if all integrals were screened out

    				for(auto f1=0, f1234=0; f1!=n1; ++f1){

						auto bf1_ind = bf1 + f1;

      				for(auto f2=0; f2!=n2; ++f2){

							auto bf2_ind = bf2 + f2;

      					for(auto f3=0; f3!=n3; ++f3){

								auto bf3_ind = bf3 + f3;

      						for(auto f4=0; f4!=n4; ++f4, ++f1234){

									auto bf4_ind = bf4 + f4;

									*(g - 1 + index_four(bf1_ind+1,bf2_ind+1,bf3_ind+1,bf4_ind+1,num_aos,num_aos,num_aos)) = ints_1234[f1234];
								}
							}
						}
					}
				}
			}
  		}
	}

	finalize();

	return;
}
