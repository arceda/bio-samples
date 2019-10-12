

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

void calculateDeltas(ublas::vector<int> individual, int i, int j, ublas::matrix<int> matrix_w, int& delta_c, int &delta_f){
    //cout<<individual<<endl;
    //cout<<matrix_w<<endl;
    //cout<<i<<" "<<j<<endl;
    //cout<<individual[i-1]<<individual[i]<<individual[j]<<individual[j+1]<<endl;
    delta_c = 0;
    delta_f = 0;
    delta_f = delta_f - matrix_w(individual[i-1], individual[i]) - matrix_w(individual[j], individual[j+1]);
    delta_f = delta_f + matrix_w(individual[i-1], individual[j]) + matrix_w(individual[i], individual[j+1]);
    if (matrix_w(individual[i-1], individual[i]) > CUTOFF)
        delta_c = delta_c + 1;
    if (matrix_w(individual[j], individual[j+1]) > CUTOFF)
        delta_c = delta_c + 1;
    if (matrix_w(individual[i-1], individual[j]) > CUTOFF)
        delta_c = delta_c - 1;
    if (matrix_w(individual[i], individual[j+1]) > CUTOFF)
        delta_c = delta_c - 1;
}

void applyMovement(ublas::vector<int>& individual, int i, int j){
    int temp = individual[i];
    individual[i] = individual[j];
    individual[j] = temp;
}

void selectMovement(ublas::matrix<int> L, int& i, int& j){
    ///////////////////////////////////////////////
    //get the posible movement with minimun delta_c
    ublas::matrix_column< ublas::matrix<int> > L_temp(L, 3);
    int min_delta_c = min_vector(L_temp);
    //cout<<L_temp[min_delta_c]<<endl;

    ublas::matrix<int> L_with_min_delta_c;
    int index_matrix = 0;
    int x = L.size1();
    for(int i = 0; i < x; i++){
        if(L(i, 3) == L_temp[min_delta_c]){
            assign_to_matrix(L_with_min_delta_c, index_matrix, 0, L(i,0));
            assign_to_matrix(L_with_min_delta_c, index_matrix, 1, L(i,1));
            assign_to_matrix(L_with_min_delta_c, index_matrix, 2, L(i,2));
            assign_to_matrix(L_with_min_delta_c, index_matrix, 3, L(i,3));
            index_matrix++;
        }
    }

    ///////////////////////////////////////////////
    //get the posible movement with maximun delta_f    
    ublas::matrix_column< ublas::matrix<int> > L_temp2(L_with_min_delta_c, 2);
    int max_delta_f = min_vector(L_temp2);
    //cout<<"mas delta f "<<L_temp2[max_delta_f]<<endl;
    ublas::matrix<int> L_with_max_delta_f;
    index_matrix = 0;
    x = L_with_min_delta_c.size1();
    for(int i = 0; i < x; i++){
        if(L_with_min_delta_c(i, 2) == L_temp2[max_delta_f]){
            assign_to_matrix(L_with_max_delta_f, index_matrix, 0, L_with_min_delta_c(i,0));
            assign_to_matrix(L_with_max_delta_f, index_matrix, 1, L_with_min_delta_c(i,1));
            assign_to_matrix(L_with_max_delta_f, index_matrix, 2, L_with_min_delta_c(i,2));
            assign_to_matrix(L_with_max_delta_f, index_matrix, 3, L_with_min_delta_c(i,3));
            index_matrix++;
        }
    }    

    i = L_with_max_delta_f(0, 0);
    j = L_with_max_delta_f(0, 1);
    //cout<<L<<endl;
    //cout<<L_with_min_delta_c<<endl;
    //cout<<L_with_max_delta_f<<" "<<i<<" "<<j<<endl;

}


ublas::vector<int> PALS(int K, ublas::vector<int> individual, ublas::matrix<int> matrix_w){
    int iterations = 0;
    while (iterations < 3000){
        ublas::matrix<int> L;
        int l_index = 0;
        int delta_c, delta_f;
        for (int i = 1; i < K; i++){
            for (int j = 0; j < K-1; j++){    
                calculateDeltas(individual, i, j, matrix_w, delta_c, delta_f);
                if (delta_c < 0 || (delta_c == 0 && delta_f > 0)){
                    assign_to_matrix(L, l_index, 0, i);
                    assign_to_matrix(L, l_index, 1, j);
                    assign_to_matrix(L, l_index, 2, delta_f);
                    assign_to_matrix(L, l_index, 3, delta_c);
                    l_index++;
                }
            }
        }
        //cout<<L<<endl;
        //break;

        
        cout<<" interation PALS: "<<iterations<<" candidates number: "<<L.size1()<<" fitness: "<<fitness(matrix_w, individual)<<" consensus: "<<consensus(matrix_w, individual)<<endl;

        iterations++;
        if (L.size1() > 0){
            //###################################################################################################
            //PALS original [1]
            int i, j;
            selectMovement(L, i, j);
            //cout<<"individual before movement:"<<individual<<endl;
            applyMovement(individual, i, j);
            //cout<<"individual after movement:"<<individual<<endl;
            //cout<<"i,j "<<i<<" "<<j<<endl;
            //break;


        }
        else
            break;

        
    }

    return individual;
}

/*
int main(){
    ublas::matrix<int> m = read_csv("../x60189_4/matrix_conservative.csv");
    int num_fragments = m.size1();
        
    ublas::vector<int> individual_t = create_shuffle_vector(num_fragments);
    
    ublas::vector<int> solution = PALS(num_fragments, individual_t, m);
    //cout<<m<<endl;    
}*/

