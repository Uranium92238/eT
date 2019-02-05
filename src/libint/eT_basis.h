
#include <libint2.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>
#include <libint2/basis.h>
#include <libint2/shell.h>


using namespace libint2;
using namespace std;

class eTBasis: public vector<BasisSet> {
   public:
      eTBasis(){}
      void add(BasisSet new_basis){
         this->push_back(new_basis);
        // cout << new_basis.nbf() << endl;
      }

      int nbf(){
         int n = 0;
         for (auto it = this->begin(); it != this->end(); ++it){
            n = n + (*it).nbf();
         }

         return n;
      }

      int max_nprim() {
         int n = 0;
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


      Shell& operator[] (int x) {
         int n = 0;
         int ind = 0;

         for (auto it = this->begin(); it != this->end(); it++){

               n += (*it).size();

               if (x < n){
                  ind = x - (n - (*it).size());
                  return (*it)[ind];
               }
         }
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
            for (int j = 0; j != atoms.size(); j++){
               if (this->operator[](i).O[0] == atoms[j].x && this->operator[](i).O[1] == atoms[j].y && this->operator[](i).O[2] == atoms[j].z){
                  result.push_back(j);
               }

            }
        }
        return result;
      }

      vector<vector<int>>  atom2shell(const vector<Atom>& atoms) {
        vector<vector<int>> result;
        int iatom = 0;
        result.resize(atoms.size());

         for (int i = 0; i != atoms.size(); i++){
            for (int j = 0; j != this->nshells(); j++){
               if (this->operator[](j).O[0] == atoms[i].x && this->operator[](j).O[1] == atoms[i].y && this->operator[](j).O[2] == atoms[i].z){
                  result[i].push_back(j);
               }

            }
        }
        return result;
      }

 };
