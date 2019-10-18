#ifndef LOCAL_SEARCH_H 
#define LOCAL_SEARCH_H

#include <iostream>     // cout, endl
#include <fstream>      // fstream
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
# IMPLEMENTACION DE PALS CON SUS MODIFICACIONES
######################################### REFERENCES ##############################################
#[1] "A New Local Search Algorithm for the DNA Fragment Assembly Problem"
#[2] "An improved problem aware local search algorithm for the DNA fragment assembly problem"
#[3] "A hybrid crow search algorithm for solving the DNA fragment assembly problem"
###################################################################################################

//compile: g++ -std=c++11 local_search.cpp utils.h  -lboost_regex
*/

void calculateDeltas(ublas::vector<int> individual, int i, int j, ublas::matrix<int> matrix_w, int& delta_c, int &delta_f);
void applyMovement(ublas::vector<int>& individual, int i, int j);
void selectMovement(ublas::matrix<int> L, int& i, int& j);
void applyMovement_PALS2many_fit(ublas::vector<int> &individual, ublas::matrix<int> L);
ublas::vector<int> PALS(int K, ublas::vector<int> individual, ublas::matrix<int> matrix_w);

#endif