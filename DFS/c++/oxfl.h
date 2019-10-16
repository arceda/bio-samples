#ifndef OXFL_H 
#define OXFL_H

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

ublas::vector<int> oxfl(ublas::vector<int> crow1, ublas::vector<int> crow2, float FL);

#endif