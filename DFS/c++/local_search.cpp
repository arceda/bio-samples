

#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <vector>
#include <string>
#include <algorithm>    // copy
#include <iterator>     // ostream_operator
#include <random>
#include <ctime> 

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

//compile: g++ -std=c++11 local_search.cpp utils.h  -lboost_regex -o local_search.out
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

//void selectMovement(ublas::matrix<int> L, int& i, int& j){
void selectMovement(ublas::matrix<int> L, int& i, int& j){
    ///////////////////////////////////////////////
    //get the posible movement with minimun delta_c
    ublas::matrix_column< ublas::matrix<int> > L_temp(L, 3);
    //std::vector<std::vector<int> > L_temp();
    int min_delta_c = min_vector(L_temp);
    //cout<<L_temp[min_delta_c]<<endl;

    ublas::matrix<int> L_with_min_delta_c;
    int index_matrix = 0;
    int x = L.size1();
    for(int i = 0; i < x; i++){
        if(L(i, 3) == L_temp[min_delta_c]){
            //assign_to_matrix(L_with_min_delta_c, index_matrix, 0, L(i,0));
            //assign_to_matrix(L_with_min_delta_c, index_matrix, 1, L(i,1));
            //assign_to_matrix(L_with_min_delta_c, index_matrix, 2, L(i,2));
            //assign_to_matrix(L_with_min_delta_c, index_matrix, 3, L(i,3));
            L_with_min_delta_c.resize(L_with_min_delta_c.size1()+1, 4);
            L_with_min_delta_c(index_matrix, 0) = L(i,0);
            L_with_min_delta_c(index_matrix, 1) = L(i,1);
            L_with_min_delta_c(index_matrix, 2) = L(i,2);
            L_with_min_delta_c(index_matrix, 3) = L(i,3);
            index_matrix++;
        }
    }

    ///////////////////////////////////////////////
    //get the posible movement with maximun delta_f    
    ublas::matrix_column< ublas::matrix<int> > L_temp2(L_with_min_delta_c, 2);
    int max_delta_f = max_vector(L_temp2);
    //cout<<"mas delta f "<<L_temp2[max_delta_f]<<endl;
    ublas::matrix<int> L_with_max_delta_f;
    index_matrix = 0;
    x = L_with_min_delta_c.size1();
    for(int i = 0; i < x; i++){
        if(L_with_min_delta_c(i, 2) == L_temp2[max_delta_f]){
            //assign_to_matrix(L_with_max_delta_f, index_matrix, 0, L_with_min_delta_c(i,0));
            //assign_to_matrix(L_with_max_delta_f, index_matrix, 1, L_with_min_delta_c(i,1));
            //assign_to_matrix(L_with_max_delta_f, index_matrix, 2, L_with_min_delta_c(i,2));
            //assign_to_matrix(L_with_max_delta_f, index_matrix, 3, L_with_min_delta_c(i,3));
            L_with_max_delta_f.resize(L_with_max_delta_f.size1()+1, 4);
            L_with_max_delta_f(index_matrix, 0) = L_with_min_delta_c(i,0);
            L_with_max_delta_f(index_matrix, 1) = L_with_min_delta_c(i,1);
            L_with_max_delta_f(index_matrix, 2) = L_with_min_delta_c(i,2);
            L_with_max_delta_f(index_matrix, 3) = L_with_min_delta_c(i,3);
            index_matrix++;
        }
    }    

    i = L_with_max_delta_f(0, 0);
    j = L_with_max_delta_f(0, 1);
    //cout<<L<<endl;
    //cout<<L_with_min_delta_c<<endl;
    //cout<<L_with_max_delta_f<<" "<<i<<" "<<j<<endl;

}

void applyMovement_PALS2many_fit(ublas::vector<int> &individual, ublas::matrix<int>& L){
    //cout<<"applyMovement_PALS2many_fit"<<endl;
    
    //sorting descending by delta_f  
    sort_by_column(L, 2);

    //cout<<"L"<<L<<endl;
    
    std::vector<int> index_used;
    int i, j, tmp;
    std::vector<int>::iterator it_i;
    std::vector<int>::iterator it_j;

    for(int posible_movement =0; posible_movement < L.size1(); posible_movement++){
        i = L(posible_movement, 0);
        j = L(posible_movement, 1);

        //cout<<"i,j "<<i<<","<<j<<"\t";

        it_i = find (index_used.begin(), index_used.end(), i);
        it_j = find (index_used.begin(), index_used.end(), j);
        //#si el movieminto no fue aplicado antes
        if(it_i == index_used.end() && it_j == index_used.end()){
            //cout<<"i, j not in vector, applied mov"<<endl;

            index_used.push_back(i);
            index_used.push_back(j);

            //swap
            tmp = individual[i];
            individual[i] = individual[j];
            individual[j] = tmp;
        }
        //else
            //cout<<"i, j are in vector"<<endl;
    }   
}

ublas::vector<int> PALS(int K, ublas::vector<int> individual, ublas::matrix<int> matrix_w){
    //ublas::vector<int> individual(individual_temp.size());
    //std::copy(individual_temp.begin(), individual_temp.end(), individual.begin());

    unsigned t0, t1;

    int iterations = 0;
    while (iterations < 90000){

        ublas::matrix<int> L;

        int l_index = 0;
        int delta_c, delta_f;

        //###################################################################################################
        //PALS original [1]
        /*
        for (int i = 1; i < K; i++){
            for (int j = 0; j < K-1; j++){    
                calculateDeltas(individual, i, j, matrix_w, delta_c, delta_f);
                if ((delta_c < 0) || (delta_c == 0 && delta_f > 0)){                
                    L.resize(L.size1()+1, 4);
                    L(l_index, 0) = i;
                    L(l_index, 1) = j;
                    L(l_index, 2) = delta_f;
                    L(l_index, 3) = delta_c;
                    l_index++;
                }
            }            
        }

        if (L.size1() > 0){
            int i, j;    
            //unsigned t0=clock();

            //selectMovement(L, i, j);
            //applyMovement(individual, i, j);

            //unsigned t1 = clock();
            //cout << "Execution Time: " << (double(t1-t0)/CLOCKS_PER_SEC) << endl;
        } else break; 
        //###################################################################################################
        */

        //###################################################################################################
        //#PALS modificado en [3]        
        for (int i = 1; i < K; i++){
            for (int j = i; j < K-1; j++){    
                calculateDeltas(individual, i, j, matrix_w, delta_c, delta_f);
                if (delta_f > 0){
                    L.resize(L.size1()+1, 4);
                    L(l_index, 0) = i;
                    L(l_index, 1) = j;
                    L(l_index, 2) = delta_f;
                    L(l_index, 3) = delta_c;
                    l_index++;
                }
            }            
        }
       
       if (L.size1() > 0){  
            //cout<<"individual before movement:"<<individual<<endl;
            //t0=clock();
            applyMovement_PALS2many_fit(individual, L);
            //t1 = clock();
            //cout << "Execution Time: " << (double(t1-t0)/CLOCKS_PER_SEC) << endl;
            //cout<<"individual after movement:"<<individual<<endl;       
            //break;  
        }else break;  
        //###################################################################################################

        //cout<<" interation PALS: "<<iterations<<" candidates number: "<<L.size1()<<endl;
        cout<<" interation PALS: "<<iterations<<" candidates number: "<<L.size1()<<endl<<L<<endl<<endl;
        //cout<<" interation PALS: "<<iterations<<" candidates number: "<<L.size1()<<" fitness: "<<fitness(matrix_w, individual)<<" consensus: "<<consensus(matrix_w, individual)<<endl;

        iterations++;      
    }

    return individual;
}


int main(){
    ublas::matrix<int> m = read_csv("../x60189_4/matrix_conservative.csv");
    int num_fragments = m.size1();
        
    ublas::vector<int> individual_t = create_shuffle_vector(num_fragments);
    
    ublas::vector<int> solution = PALS(num_fragments, individual_t, m);
    cout<<"fit :"<<fitness(m, solution)<<endl;    
    cout<<"contigs :"<<consensus(m, solution)<<endl;    


    /*
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////PRUEBAS///////////////////////////////////////////////////
    cout<<"TESTING... :"<<endl<<endl;    
   
    int f_tmp, c_tmp;
    std::vector<int> fitness_vec;
    std::vector<int> contigs_vec;
    for (int i = 0; i < 30; i++){
        ublas::vector<int> individual = create_shuffle_vector(num_fragments);    
        ublas::vector<int> sol = PALS(num_fragments, individual, m);
        f_tmp = fitness(m, sol);
        c_tmp = consensus(m, sol);
        cout<<"fit :"<<f_tmp<<"\t"<<"contigs :"<<c_tmp<<endl;  
        fitness_vec.push_back(f_tmp);
        contigs_vec.push_back(c_tmp);
        
    }
    int best_fitness = *std::max_element(fitness_vec.begin(), fitness_vec.end());
    int worst_fitness = *std::min_element(fitness_vec.begin(), fitness_vec.end());
    float mean_fitness = (std::accumulate(fitness_vec.begin(), fitness_vec.end(), 0.0))/fitness_vec.size();

    int best_contigs = *std::min_element(contigs_vec.begin(), contigs_vec.end());
    int worst_contigs = *std::max_element(contigs_vec.begin(), contigs_vec.end());
    float mean_contigs = (std::accumulate(contigs_vec.begin(), contigs_vec.end(), 0.0))/contigs_vec.size();

    cout<<endl;
    cout<<"BEST"<<"\t\t"<<"MEAN"<<"\t\t"<<"WORST"<<endl;
    cout<<best_fitness<<"\t\t"<<mean_fitness<<"\t\t"<<worst_fitness<<endl;
    cout<<best_contigs<<"\t\t"<<mean_contigs<<"\t\t"<<worst_contigs<<endl;    
    ///////////////////////////////////////////////////////////////////////////////////////////////  
    */
}

