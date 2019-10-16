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

        //cout<<"cecee "<<ublas::row(crows, i)<<endl;
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
    int random_crow_index;
    float r, r_ls;
    ublas::vector<int> current_crow;
    ublas::vector<int> current_memory;
    while(iter < 300){
        cout<<"ITERATION "<<iter<<endl;
        for (int i = 0; i < N; i++){
            cout<<"crow "<<i<<endl;
            current_crow = row(crows, i);
            current_memory = row(memory, i);
            //cout<<"current_crow "<<current_crow<<endl;
            //cout<<random_number(0, N-1)<<endl;
            r = float(random_float());            
            if (r >= AP){
                random_crow_index = random_number(0, N-1);
                //OXFL oeprator          
                
                r_ls = float(random_float());
                if (r_ls >= P_LS){
                    cout<<"local search ..."<<endl;
                    ublas::row(crows, i) = PALS(num_fragments, current_crow, matrix_w);    
                    current_crow =   ublas::row(crows, i); 
                }
            }else{
                //move to random position     
                ublas::row(crows, i) = create_shuffle_vector(num_fragments); 
                current_crow =   ublas::row(crows, i); 
            }

            if (fitness(matrix_w, current_crow) > fitness(matrix_w, current_memory)){
                ublas::vector<int> temp(current_crow.size());
                std::copy(current_crow.begin(), current_crow.end(), temp.begin());

                ublas::row(memory, i) = temp;                 
            }
          
        }
        iter++;
    }
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