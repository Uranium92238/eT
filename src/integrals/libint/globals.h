/*
/
/ 	Globals
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
/ 	Contains the only global variables allowed in the integral calculations,
/ 	namely the basis set, atomic information, and integral engines.
/
*/
extern libint2::BasisSet basis; 					// The basis set used throughout
extern libint2::Engine electronic_repulsion; // The electronic repulsion engine
extern libint2::Engine kinetic; 					// The kinetic energy engine
extern libint2::Engine nuclear; 					// The nuclear attraction engine
extern libint2::Engine overlap; 					// The AO overlap engine
extern vector<libint2::Atom> atoms; 			// Atoms vector
