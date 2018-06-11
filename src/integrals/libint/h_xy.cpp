/*
/
/ 	One-electron integral routines (for h_xy = h_αβ)
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
*/
#include "utils.h"
#include "h_xy.h"

#include <libint2.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace libint2;
using namespace std;

void get_n_aos(int *n_ao){
//
	initialize();
//
	string xyzfilename = "Water.xyz"; // see XYZ format description at http://en.wikipedia.org/wiki/XYZ_file_format
	ifstream input_file(xyzfilename);
	vector<Atom> atoms = read_dotxyz(input_file);
//
	BasisSet obs("cc-pVDZ", atoms);
//
	*n_ao = 0;
//
	for(auto s1=0; s1!=obs.size(); ++s1){
//
		*n_ao = *n_ao + obs[s1].size(); // Add number of basis functions in the given shell
//
	}
//
	finalize();
//
	return;
//
}

void get_ao_xy_kinetic(double *h){
//
// ** Initialize libint calculator **
//
	initialize();
//
// ** Read molecular geometry and make atoms array **
//
	string xyzfilename = "Water.xyz"; // see XYZ format description at http://en.wikipedia.org/wiki/XYZ_file_format
	ifstream input_file(xyzfilename);
	vector<Atom> atoms = read_dotxyz(input_file);
//
// ** Set basis set to use for the atoms **
//
	BasisSet obs("cc-pVDZ", atoms);
//
// ** Calculate number of basis functions from the basis object
//
	int num_aos = 0;
//
	for(auto s1=0; s1!=obs.size(); ++s1){
//
		num_aos = num_aos + obs[s1].size(); // Add number of basis functions in the given shell
//
	}
//
	cout << "There are " << num_aos << " number of basis functions for H2O/cc-pVDZ" << endl;
//
// ** Compute the kinetic energy part of one-electron integrals **
//
	Engine k_engine(Operator::kinetic,  // will compute kinetic ints overlap
                	obs.max_nprim(),     // max # of primitives in shells this engine will accept
                	obs.max_l()          // max angular momentum of shells this engine will accept
               	);
//
	auto shell2bf = obs.shell2bf(); // maps shell index to basis function index
                                   // shell2bf[0] = index of the first basis function in shell 0
                                   // shell2bf[1] = index of the first basis function in shell 1
                                   // ...
//
	const auto& buf_vec = k_engine.results(); // will point to computed shell sets
//
	int counter = 0; /// must be smarter with the calculation of position... this is incorrect...
//
	for(auto s1=0; s1!=obs.size(); ++s1) {
  		for(auto s2=0; s2!=obs.size(); ++s2) {
//
    		cout << "compute shell set {" << s1 << "," << s2 << "} ... ";
    		k_engine.compute(obs[s1], obs[s2]);
    		cout << "done" << endl;
    		auto ints_shellset = buf_vec[0];  // location of the computed integrals
    		if (ints_shellset == nullptr)
      		continue;  // nullptr returned if the entire shell-set was screened out

    		auto bf1 = shell2bf[s1];  // first basis function in first shell
    		auto n1 = obs[s1].size(); // number of basis functions in first shell
    		auto bf2 = shell2bf[s2];  // first basis function in second shell
    		auto n2 = obs[s2].size(); // number of basis functions in second shell

    		// integrals are packed into ints_shellset in row-major (C) form
    		// this iterates over integrals in this order
    		for(auto f1=0; f1!=n1; ++f1){
      		for(auto f2=0; f2!=n2; ++f2){
					*(h + counter) = *(h + counter) + ints_shellset[f1*n2+f2];
					counter = counter + 1;
    //    			cout << "  " << bf1+f1 << " " << bf2+f2 << " " << ints_shellset[f1*n2+f2] << endl;
				}
			}
  		}
	}
//
//	** Compute the nuclear attraction energy part of one-electron integrals **
//
	Engine n_engine(Operator::nuclear,  // will compute kinetic ints overlap
                	obs.max_nprim(),     // max # of primitives in shells this engine will accept
                	obs.max_l()          // max angular momentum of shells this engine will accept
               	);
//
// Make point charges (R_I, Z_I)
//
	n_engine.set_params(make_point_charges(atoms));  // convert `atoms` to point charges

//
	const auto& buf_vec_n = n_engine.results(); // will point to computed shell sets
                                          	  // const auto& is very important!
//
	counter = 0;
//
	for(auto s1=0; s1!=obs.size(); ++s1) {
  		for(auto s2=0; s2!=obs.size(); ++s2) {

    		cout << "compute shell set {" << s1 << "," << s2 << "} ... ";
    		n_engine.compute(obs[s1], obs[s2]);
    		cout << "done" << endl;
    		auto ints_shellset = buf_vec_n[0];  // location of the computed integrals
    		if (ints_shellset == nullptr)
      		continue;  // nullptr returned if the entire shell-set was screened out

    		auto bf1 = shell2bf[s1];  // first basis function in first shell
    		auto n1 = obs[s1].size(); // number of basis functions in first shell
    		auto bf2 = shell2bf[s2];  // first basis function in second shell
    		auto n2 = obs[s2].size(); // number of basis functions in second shell

    		// integrals are packed into ints_shellset in row-major (C) form
    		// this iterates over integrals in this order
    		for(auto f1=0; f1!=n1; ++f1){
      		for(auto f2=0; f2!=n2; ++f2){
        			cout << "  " << bf1+f1 << " " << bf2+f2 << " " << ints_shellset[f1*n2+f2] << endl;
			//		*(h + counter) = *(h + counter) + ints_shellset[f1*n2+f2];
					counter = counter + 1;
				}
			}
//
  		}
	}
//
	cout << index_two(1, 1, 2) << " should be zero!" << endl;
	cout << index_two(2, 1, 2) << " should be one!" << endl;
//
// ** Finalize libint calculator **
//
	finalize();
//
	return;
}
