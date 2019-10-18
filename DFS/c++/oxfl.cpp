#include <iostream>     // cout, endl
#include <vector>
#include <string>
#include <algorithm>    // copy
#include <iterator>     // ostream_operator
#include <random>

#include <boost/tokenizer.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>


#include "utils.h"

using namespace std;
using namespace boost; //range
using namespace boost::numeric;//ublas.matrix, ublas.vector
using namespace boost::numeric::ublas;//matrix

/*
# IMPLEMENTACION DEL OPERADOR OXFL
######################################### REFERENCES ##############################################
#[3] "A hybrid crow search algorithm for solving the DNA fragment assembly problem"
###################################################################################################

//compile: g++ -std=c++11 oxfl.cpp utils.h  -lboost_regex
*/

ublas::vector<int> operator_oxfl(ublas::vector<int> crow, ublas::vector<int> victim, float FL){
    
    int n = crow.size();

    //C1 y C2 van de 0 a num_fragments-1
    //int C1 = int(random_number(0, n-1))
    //int C2 = int((C1+(n-1)*FL) % (n-1))
    int C1 = 11;
    int C2 = 1; 
    int temp;

    ublas::vector<int> sol(n);
    std::vector<int> X(n);
    std::vector<int> M(n);

    if (C1 < C2){
        //X = crow;
        std::copy(crow.begin(), crow.end(), X.begin());
        //M = victim;   
        std::copy(victim.begin(), victim.end(), M.begin());
    }else{
        //X = victim;
        std::copy(victim.begin(), victim.end(), X.begin());
        //M = crow;
        std::copy(crow.begin(), crow.end(), M.begin());
        temp = C1;
        C1 = C2;
        C2 = temp;
    }

    //sol[C1:C2] = X[C1:C2]
    for (int i = C1; i < C2; i++){
        sol[i] = X[i];
    }

    //delete the range C1:c2 from M
    for (int i = C1; i < C2; i++){
        std::vector<int>::iterator index = std::find(M.begin(), M.end(), X[i]);
        M.erase(index);
    }

    //sol[0:C1] = M[0:C1] # insert the firsts to C1 to new solution
    for (int i = 0; i < C1; i++){
        sol[i] = M[i];
        M.erase(M.begin() + i);
    }

    //insert he ñast elents of M
    for (int i = 0; i < M.size(); i++){
        sol[C2 + i] = M[i];
    }

    //cout<<"C1:"<<C1<<endl;
    //cout<<"C2:"<<C2<<endl;
    //cout<<"crow:"<<crow<<endl;    
    //cout<<"victim:"<<victim<<endl;
    //cout<<"sol:"<<sol<<endl;
    
    //for(int i = 0; i < M.size(); i++)
    //    cout<<"M_"<<i<<":"<<M[i]<<endl;

    return sol;   
}

/*
int main(){

    std::vector<int> tmp{9,11,5,8,2,13,14,1,4,3,7,10,6,12};
    ublas::vector<int> crow(tmp.size());
    std::copy(tmp.begin(), tmp.end(), crow.begin());

    std::vector<int> tmp2{4,10,7,2,13,8,11,5,14,12,6,1,3,9};
    ublas::vector<int> victim(tmp2.size());
    std::copy(tmp2.begin(), tmp2.end(), victim.begin());
   
    float FL = 0.75;
    
    ublas::vector<int> sol = operator_oxfl(crow, victim, FL);    
    
  
}*/