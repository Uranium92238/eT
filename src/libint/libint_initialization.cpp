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

    Libint initialization of engines and basis set
    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018

*/
#include "libint_initialization.h"
#include "omp_control.h"

#include "eT_basis.h"

#include <libint2.hpp>

using namespace libint2;
using namespace std;

// Global variables declared in this file

// Basis set and atoms
eTBasis basis;
vector<Atom> atoms;

/*
   Vectors of Libint engines.
   Length of vectors is equal to number of threads
*/
vector<Engine> electronic_repulsion(omp_get_max_threads());
vector<Engine> electronic_repulsion_1der(omp_get_max_threads());
vector<Engine> kinetic(omp_get_max_threads());
vector<Engine> kinetic_1der(omp_get_max_threads());
vector<Engine> nuclear(omp_get_max_threads());
vector<Engine> nuclear_1der(omp_get_max_threads());
vector<Engine> overlap(omp_get_max_threads());
vector<Engine> overlap_1der(omp_get_max_threads());
vector<Engine> dipole(omp_get_max_threads());
vector<Engine> quadrupole(omp_get_max_threads());
vector<Engine> coulomb_external_charges(omp_get_max_threads());
vector<Engine> coulomb_external_unit_charges(omp_get_max_threads());


void export_geometry_and_basis_to_libint(const int nAtoms,
                                         const int *atomicNumbers,
                                         const double *atomicCoordinates,
                                         const char *basisSets,
                                         const int maxLength,
                                         const int *cartesians){

    /*
    Export geometry and basis set
    Written by Rolf H. Myhre, Mar. 2020

    Exports the geometry and basis sets used from the Fortran side to Libint

    n_atoms is an int with number of atoms

    atomic_numbers is an int array of dimension n_atoms passed by reference

    atomic_cordinates is a double array with dimension 3*n_atoms with atomic coordinates
    passed by reference

    basisSets is an array of c strings passed by reference

    maxLength is an integer indicating the storage length of the strings in basisSet

    cartesians is an array of ints passed by reference
    */

    //Allocate a temporary atom object
    vector<Atom> temporaryAtoms(nAtoms, Atom{0, 0.0, 0.0, 0.0});

    //Allocate vectors to keep track of atoms and basis sets
    vector<string> uniqueBasis;
    vector<int> uniqueCartesians;
    vector<vector<int>> atomsWithBasis;

    //Some extra variables
    bool newBasis;
    int nBasis, currentBasis;

    //Construct std::string using the reference to the first basis set
    //and add to list of unique basis sets
    string basisString = string(&basisSets[0]);
    uniqueBasis.push_back(basisString);

    //Initialize vectors with cartesians and list of atom with each basis
    uniqueCartesians.push_back(cartesians[0]);
    atomsWithBasis.push_back({});
    nBasis = 1;

    //Loop over the atoms
    for (int i = 0; i < nAtoms; i++)
    {
        //set atomic number and coordinates in temporaryAtoms coordinates
        temporaryAtoms[i].atomic_number = atomicNumbers[i];
        temporaryAtoms[i].x = atomicCoordinates[3*i];
        temporaryAtoms[i].y = atomicCoordinates[3*i+1];
        temporaryAtoms[i].z = atomicCoordinates[3*i+2];

        //Construct std::string by referencing basisSets offset by maxLength*i
        basisString = string(&basisSets[i*maxLength]);

        //Check if there's a new basis
        newBasis = true;
        for (int j = 0; j < nBasis; j++){
            if(uniqueBasis[j] == basisString &&
               uniqueCartesians[j] == cartesians[i]){

                newBasis = false;
                currentBasis = j;
                break;
            }
        }

        //If we found a new basis, add it to the list,
        //else just add the index to the list for the current basis
        if(newBasis){

            nBasis++;
            uniqueBasis.push_back(basisString);
            uniqueCartesians.push_back(cartesians[i]);
            atomsWithBasis.push_back({i});
        }
        else{
            atomsWithBasis[currentBasis].push_back(i);
        }
    }

    //Set atoms in Libint
    atoms = temporaryAtoms;

    //Reset basis in libint
    eTBasis tempBasis;
    basis = tempBasis;

    //Set up the new temporary basis
    BasisSet temporary;

    for(int i = 0; i < nBasis; i++){

        //clear temporary atoms
        temporaryAtoms.clear();

        //Add atoms with uniqueBasis i to temporaryAtoms
        for (int const& atom: atomsWithBasis[i]){
            temporaryAtoms.push_back(atoms[atom]);
        }

        //Create temporary basis set for uniqueBasis i
        temporary = BasisSet(uniqueBasis[i], temporaryAtoms, true);

        //Set pure or cartesian for uniqueBasis i
        for(auto& shell: temporary) {
            for(auto& contraction: shell.contr) {

                if (contraction.l >= 2 && uniqueCartesians[i] != 0){
                    contraction.pure = false;
                }
            }
        }

        //Add temporary to Libint basis
        basis.add(temporary);
    }
}



void initialize_libint(){

    initialize();

}

void finalize_libint(){

    finalize();

}

void initialize_eri(const double eri_precision){

    electronic_repulsion[0] = Engine(Operator::coulomb, basis.max_nprim(), basis.max_l());
    electronic_repulsion_1der[0] = Engine(Operator::coulomb, basis.max_nprim(), basis.max_l(), 1);

    electronic_repulsion[0].set_precision(eri_precision);
    electronic_repulsion_1der[0].set_precision(eri_precision);

    for (int i = 1; i < omp_get_max_threads(); i++){

        electronic_repulsion[i] = electronic_repulsion[0];
        electronic_repulsion_1der[i] = electronic_repulsion_1der[0];

    }

}

void set_eri_precision(const double eri_precision){

    for (int i = 0; i < omp_get_max_threads(); i++){

        electronic_repulsion[i].set_precision(eri_precision);
        electronic_repulsion_1der[i].set_precision(eri_precision);

    }

}

void initialize_kinetic(){


    kinetic[0] = Engine(Operator::kinetic, basis.max_nprim(), basis.max_l());
    kinetic_1der[0] = Engine(Operator::kinetic, basis.max_nprim(), basis.max_l(),1);

    for (int i = 1; i < omp_get_max_threads(); i++){

        kinetic[i] = kinetic[0];
        kinetic_1der[i] = kinetic_1der[0];

    }

}

void initialize_nuclear(){

    vector<pair<double, array<double, 3>>> q;
    for (const auto& atom : atoms) {
      q.push_back({static_cast<double>(atom.atomic_number),
                   {{atom.x, atom.y, atom.z}}});
    }

    nuclear[0] = Engine(Operator::nuclear, basis.max_nprim(), basis.max_l());
    nuclear_1der[0] = Engine(Operator::nuclear, basis.max_nprim(), basis.max_l(), 1);

    nuclear[0].set_params(q);
    nuclear_1der[0].set_params(q);

    for (int i = 1; i < omp_get_max_threads(); i++){

        nuclear[i] = nuclear[0];
        nuclear_1der[i] = nuclear_1der[0];

    }

}

void initialize_overlap(){

    overlap[0] = Engine(Operator::overlap, basis.max_nprim(), basis.max_l());
    overlap_1der[0] = Engine(Operator::overlap, basis.max_nprim(), basis.max_l(),1);

    for (int i = 1; i < omp_get_max_threads(); i++){

        overlap[i] = overlap[0];
        overlap_1der[i] = overlap_1der[0];

    }

}

void initialize_dipole(){

    dipole[0] = Engine(Operator::emultipole1, basis.max_nprim(), basis.max_l()); // Note that expansion point = (0,0,0) by default

    for (int i = 1; i < omp_get_max_threads(); i++){
        dipole[i] = dipole[0];
    }

}

void initialize_quadrupole(){

    quadrupole[0] = Engine(Operator::emultipole2, basis.max_nprim(), basis.max_l()); // Note that expansion point = (0,0,0) by default

    for (int i = 0; i < omp_get_max_threads(); i++){
        quadrupole[i] = quadrupole[0];
    }

}

void initialize_coulomb_external_charges(const double *charges, const double *coordinates,
                                         const int n_points){

    coulomb_external_charges[0] = Engine(Operator::nuclear, basis.max_nprim(), basis.max_l());

    vector<pair<double, array<double, 3>>> points;

    for (int i = 0; i < n_points; i++) {

      int offset = (i-1)*3 + 3;
      points.push_back({(charges[i]),
                   {{*(coordinates+offset),*(coordinates+offset+1),
                    *(coordinates+offset+2)}}});

    }

    coulomb_external_charges[0].set_params(points);

    for (int i = 1; i < omp_get_max_threads(); i++){

        coulomb_external_charges[i] = coulomb_external_charges[0];

    }

}

void initialize_coulomb_external_unit_charges(const double *coordinates, const int n_points){

    coulomb_external_unit_charges[0] = Engine(Operator::nuclear, basis.max_nprim(), basis.max_l());

    vector<pair<double, array<double, 3>>> points;

    for (int i = 0; i < n_points; i++) {

      int offset = (i-1)*3 + 3;
      points.push_back({static_cast<double>(1.0e0),
                   {{*(coordinates+offset),*(coordinates+offset+1),
                    *(coordinates+offset+2)}}});
    }

    coulomb_external_unit_charges[0].set_params(points);

    for (int i = 1; i < omp_get_max_threads(); i++){

        coulomb_external_unit_charges[i] = coulomb_external_unit_charges[0];

    }

}
