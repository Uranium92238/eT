//
//
//  eT - a coupled cluster program
//  Copyright (C) 2016-2022 the authors of eT
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
 
    eT basis
    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018

    Enables different basis sets for different atoms

*/
#include <libint2.hpp>
#include <libint2/basis.h>
#include <libint2/shell.h>

using namespace libint2;
using namespace std;

class eTBasis: public vector<BasisSet> {
   public:
      eTBasis(){}
      void add(BasisSet new_basis){
         this->push_back(new_basis);
      }

      int nbf(){
         int n = 0;
         for (auto it = this->begin(); it != this->end(); ++it){
            n = n + (*it).nbf();
         }

         return n;
      }

      int max_nprim() {
         std::size_t n = 0;
         for (auto it = this->begin(); it != this->end(); ++it){
            if ((*it).max_nprim() > n){
               n = (*it).max_nprim();
            }
            
         }
        return n;
      }

      int max_l() {
         int n = 0;
         for (auto it = this->begin(); it != this->end(); ++it){
            if ((*it).max_l() > n){
               n = (*it).max_l();
            }
            
         }
        return n;
      }

      int nshells() {
         int n = 0;
         for (auto it = this->begin(); it != this->end(); it++){
               n += (*it).size(); 
         }
        return n;
      }


      const Shell& operator[] (int x) {
         int n = 0;
         int ind = 0;

         for (auto it = this->begin(); it != this->end(); it++){

               n += (*it).size();

               if (x < n){
                  ind = x - (n - (*it).size());
                  return (*it)[ind];
               }
         }
         throw std::invalid_argument("x too big");
      }

      vector<int> shell2bf() {
         vector<int> result;

         result.reserve(this->nshells());

         int n = 0;
          for (auto i = 0; i != this->nshells(); i++){

             result.push_back(n);
             n += this->operator[](i).size();
          }

         return result;
      }

      vector<int> shell2atom(const vector<Atom>& atoms) {
        vector<int> result;
        result.reserve(nshells());

         for (int i = 0; i != this->nshells(); i++){
            for (std::size_t j = 0; j != atoms.size(); j++){
               if (this->operator[](i).O[0] == atoms[j].x && this->operator[](i).O[1] == atoms[j].y && this->operator[](i).O[2] == atoms[j].z){
                  result.push_back(j);
               }

            }
        }
        return result;
      }

      vector<vector<int>>  atom2shell(const vector<Atom>& atoms) {
        vector<vector<int>> result;

        result.resize(atoms.size());

         for (std::size_t i = 0; i != atoms.size(); i++){
            for (int j = 0; j != this->nshells(); j++){
               if (this->operator[](j).O[0] == atoms[i].x && this->operator[](j).O[1] == atoms[i].y && this->operator[](j).O[2] == atoms[i].z){
                  result[i].push_back(j);
               }

            }
        }
        return result;
      }

 };
