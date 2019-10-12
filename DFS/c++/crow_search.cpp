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
#include "local_search.cpp"

using namespace std;
using namespace boost; //range
using namespace boost::numeric;//ublas.matrix, ublas.vector
using namespace boost::numeric::ublas;//matrix

//compile: g++ -std=c++11 crow_search.cpp  local_search.h utils.h  -lboost_regex



ublas::matrix<int> init_population(int N, int num_fragments){
    ublas::matrix<int> crows(N, num_fragments);
    //cout<<"crows "<<crows<<endl;
    for(int i = 0; i < N; i++){
        ublas::vector<int> individual = create_shuffle_vector(num_fragments);
        //cout<<"individual "<<individual<<endl;
        ublas::row(crows, i) = individual;        
    }
    //cout<<"crows "<<crows<<endl;

    return crows;
}

/*def init_population():
    for i in range(N):
        crow = np.arange(num_fragments) #individual = [0, 1, 2, ...] each index is a fragment
        np.random.shuffle(crow) #shuffle the fragment, this a ramdon solution
        crows[i] = crow
    return crows
*/
void crow_search(ublas::matrix<int> matrix_w, int ITERATIONS, int N, float AP, float FL, float P_LS){
    int num_fragments = matrix_w.size1();

    ublas::matrix<int> crows = init_population(N, num_fragments);
    ublas::matrix<int> memory(crows);

    int iter = 0;
    int r, r_ls;
    while(iter < 300){
        for (int i = 0; i < N; i++){
            cout<<"CROW "<<i<<endl;
            //cout<<random_number(0, N-1)<<endl;
            r = float(random_number(0, N-1));
            if (r >= AP){

            }else{

            }
        }
        iter++;
    }

    //cout<<"crows "<<crows<<endl;
    //cout<<"memory "<<memory<<endl;

/*
    iter=0
while iter < ITERATIONS:
    print("ITERATION: ", iter)
    for i in range(N):
        print('CROW ', i)
        random_crow = random.randint(0, N-1) #chose a random crow
        r = random.random()
        if r >= AP:
            #print("the crow look up", i)
            #print("perform oxfl operator")
            crows[i] = OXFL(crows[i], crows[random_crow], FL)

            #################     local search     ###################
            r_ls = random.random()
            if r_ls >= P_LS:
                print("local search...")
                individual = crows[i].copy()
                individual = np.squeeze(np.asarray(individual))
                crows[i] = PALS(num_fragments, individual, matrix)                

        else:
            #print("the crow move to ramdon position", i)
            #the crow go to a random position
            #print('the crow go to random position', i, crows[i])
            np.random.shuffle(crows[i])
            #print('the crow went to random position', i, crows[i])
            #print("memory[i]: ",i,  memory[i])
            #print("crows[i]: ",i,  crows[i])

        if fitness(matrix, crows[i]) > fitness(matrix, memory[i]):
            #print("the new position is better, updating memory")
            memory[i] = crows[i].copy()
            #print("memory[i]: ",i,  memory[i])
            #print("crows[i]: ",i,  crows[i])
        
    iter += 1
*/
}


int main(){
    int ITERATIONS = 200;
    int N = 32;
    float AP = 0.02;
    float FL = 0.75;
    float P_LS = 0.49;

    ublas::matrix<int> m = read_csv("../x60189_4/matrix_conservative.csv");
    int num_fragments = m.size1();
    
    crow_search(m, ITERATIONS, N, AP, FL, P_LS);

}