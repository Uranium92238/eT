/*
/
/ 	One-electron integral routines (for s_xy = s_αβ)
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
*/
#include <iostream>

using namespace std;

#include <fstream>
#include <string>
#include <vector>
#include "utils.h"
#include "s_xy.h"

#include <libint2.hpp>

#include "globals.h"

using namespace libint2;

void get_n_shells(int *ns){

	//initialize();

	*ns = basis.size();

	//finalize();

	return;

}

void get_ao_s_xy(double *s){

	//initialize();

	int num_aos = 0;
	num_aos = basis.nbf();

	auto shell2bf = basis.shell2bf(); // maps shell index to basis function index
                                   // shell2bf[0] = index of the first basis function in shell 0
                                   // shell2bf[1] = index of the first basis function in shell 1
                                   // ...

	const auto& buf_vec = overlap.results(); // will point to computed shell sets

	for(auto s1=0; s1!=basis.nshells(); ++s1) {
  		for(auto s2=0; s2!=basis.nshells(); ++s2) {

    		overlap.compute(basis[s1], basis[s2]);
    		auto ints_shellset = buf_vec[0];  // location of the computed integrals
    		if (ints_shellset == nullptr)
      		continue;  // nullptr returned if the entire shell-set was screened out

    		auto bf1 = shell2bf[s1];  // first basis function in first shell
    		auto n1 = basis[s1].size(); // number of basis functions in first shell
    		auto bf2 = shell2bf[s2];  // first basis function in second shell
    		auto n2 = basis[s2].size(); // number of basis functions in second shell

    		// integrals are packed into ints_shellset in row-major (C) form
    		// this iterates over integrals in this order
    		for(auto f1=0; f1!=n1; ++f1){
      		for(auto f2=0; f2!=n2; ++f2){
					*(s - 1 + index_two(bf1+1+f1, bf2+1+f2, num_aos)) = ints_shellset[f1*n2+f2];
				}
			}
  		}
	}

	//finalize();

	return;
}
