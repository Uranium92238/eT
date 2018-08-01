
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
         cout << new_basis.nbf() << endl;
      }

      long nbf(){
         long n = 0;
         for (auto it = this->begin(); it != this->end(); ++it){
            n = n + (*it).nbf();
         }

         return n;
      }

      long max_nprim() {
         long n = 0;
         for (auto it = this->begin(); it != this->end(); ++it){
            if ((*it).max_nprim() > n){
               n = (*it).max_nprim();
            }
            
         }
        return n;
      }

      long max_l() {
         long n = 0;
         for (auto it = this->begin(); it != this->end(); ++it){
            if ((*it).max_l() > n){
               n = (*it).max_l();
            }
            
         }
        return n;
      }

      long nshells() {
         long n = 0;
         for (auto it = this->begin(); it != this->end(); it++){
               n += (*it).size(); 
         }
        return n;
      }


      Shell& operator[] (int x) {
         long n = 0;
         long ind = 0;

         for (auto it = this->begin(); it != this->end(); it++){

               n += (*it).size();

               if (x < n){
                  ind = x - (n - (*it).size());
                  return (*it)[ind];
               }
         }
      }

      vector<long> shell2bf() {
         vector<long> result;

         result.reserve(this->nshells());

         long n = 0;
          for (auto i = 0; i != this->nshells(); i++){

             result.push_back(n);
             n += this->operator[](i).size();
          }

         return result;
      }

      vector<long> shell2atom(const vector<Atom>& atoms) {
        vector<long> result;
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

      vector<vector<long>>  atom2shell(const vector<Atom>& atoms) {
        vector<vector<long>> result;
        long iatom = 0;
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