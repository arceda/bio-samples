#include <iostream>     // cout, endl
#include <vector>
#include <string>
#include <algorithm>    // copy
#include <iterator>     // ostream_operator
#include <random>

#include <boost/tokenizer.hpp>
#include <boost/range/algorithm.hpp>
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

ublas::vector<int> operator_oxfl(ublas::vector<int> crow1, ublas::vector<int> crow2, float FL){
    crow1[0] = 100;

    cout<<"crow1"<<crow1<<endl;  
    return crow2;
}

int main(){

    std::vector<int> tmp{9,11,5,8,2,13,14,1,4,3,7,10,6,12};
    ublas::vector<int> crow(tmp.size());
    std::copy(tmp.begin(), tmp.end(), crow.begin());

    std::vector<int> tmp2{4,10,7,2,13,8,11,5,14,12,6,1,3,9};
    ublas::vector<int> victim(tmp2.size());
    std::copy(tmp2.begin(), tmp2.end(), victim.begin());

    //ublas::vector<int> crow(14);
    //crow[0] = 9;crow[1] = 11;crow[2] = 5;crow[3] = 8;crow[4] = 2;

    //ublas::vector<int> victim(14);
    //crow[0] = 9;crow[1] = 11;crow[2] = 5;crow[3] = 8;crow[4] = 2;

    int C1 = 2;
    int C2 = 11;
    float FL = 0.75;

    cout<<crow<<endl;  
    ublas::vector<int> sol = operator_oxfl(crow, victim, FL);    
    cout<<crow<<endl;    
  
}