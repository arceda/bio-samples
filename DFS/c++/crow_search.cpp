#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <vector>
#include <string>
#include <algorithm>    // copy
#include <iterator>     // ostream_operator
#include <random>
#include <ctime> 
#include <fstream>

#include <boost/tokenizer.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>


#include "utils.h"
#include "local_search.cpp"
#include "oxfl.cpp"

using namespace std;
using namespace boost; //range
using namespace boost::numeric;//ublas.matrix, ublas.vector
using namespace boost::numeric::ublas;//matrix

//compile: g++ -std=c++11 oxfl.h crow_search.cpp  local_search.h utils.h  -lboost_regex -o crow_search.out



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

void best_fitness(ublas::matrix<int> &matrix_w, ublas::matrix<int> &memory, int &fit, int &index){
    fit = 0;
    int best_index;
    int x = memory.size1();
    for (int i = 0; i < x; i++){
        int tmp_fit = fitness(matrix_w, ublas::row(memory, i));
        //int tmp_contigs = consensus(matrix_w, ublas::row(memory, i));
               
        if(tmp_fit > fit){
            fit = tmp_fit;
            index = i;            
        }            
    }
}

/*def init_population():
    for i in range(N):
        crow = np.arange(num_fragments) #individual = [0, 1, 2, ...] each index is a fragment
        np.random.shuffle(crow) #shuffle the fragment, this a ramdon solution
        crows[i] = crow
    return crows
*/
ublas::vector<int> crow_search(ublas::matrix<int> matrix_w, int ITERATIONS, int N, float AP, float FL, float P_LS, bool save=false){
    int num_fragments = matrix_w.size1();

    ublas::matrix<int> crows = init_population(N, num_fragments);
    ublas::matrix<int> memory(crows);
    
    int iter = 0;
    int random_crow_index;
    float r, r_ls;
    ublas::vector<int> current_crow;
    ublas::vector<int> current_victim;
    ublas::vector<int> current_memory;

    ofstream myfile;
    if(save) myfile.open ("_x60189_7_fitness_by_iteration.txt");
    
    

    while(iter < ITERATIONS){
        cout<<"ITERATION "<<iter<<endl;
        for (int i = 0; i < N; i++){
            //cout<<"crow "<<i<<endl;
            current_crow = ublas::row(crows, i);
            current_memory = ublas::row(memory, i);
            //cout<<"current_crow "<<current_crow<<endl;
            //cout<<random_number(0, N-1)<<endl;
            r = float(random_float());            
            if (r >= AP){
                //cout<<"oxfl ..."<<endl;
                random_crow_index = random_number(0, N-1);
                current_victim = ublas::row(crows, random_crow_index);
                //OXFL oeprator   
                ublas::row(crows, i) = operator_oxfl(current_crow, current_victim, FL); 
                current_crow =   ublas::row(crows, i);       
                
                r_ls = float(random_float());
                if (r_ls >= P_LS){
                    //cout<<"local search ..."<<endl;
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
        ////////////////////////////////////////////////////////////////////////        
        if(save) {
            int best_index, fit;
            best_fitness(matrix_w, memory, fit, best_index);
            myfile <<iter<<", "<<fit<< ", " << ublas::row(memory, best_index)<<"\n";
            cout <<iter<<", "<<fit<< ", " << ublas::row(memory, best_index)<<"\n";
        }        
        ////////////////////////////////////////////////////////////////////////


        iter++;
    }

    myfile.close();    
    /////////////////////////////////////////////////////////////////////////////////
    int best_index, fit;
    best_fitness(matrix_w, memory, fit, best_index); 
    return ublas::row(memory, best_index);
}


int main(){
    int ITERATIONS = 500;
    int N = 32;
    float AP = 0.02;
    float FL = 0.75;
    float P_LS = 0.49;

    ublas::matrix<int> m = read_csv("../x60189_7/matrix_conservative.csv");
    int num_fragments = m.size1();    
    crow_search(m, ITERATIONS, N, AP, FL, P_LS, true);

    //////////////////////////////////////////////////////////////////////////////////////////////
    /*cout<<"TESTING... :"<<endl<<endl;    
   
    int f_tmp, c_tmp;
    std::vector<int> fitness_vec;
    std::vector<int> contigs_vec;
    for (int i = 0; i < 30; i++){
        ublas::vector<int> sol = crow_search(m, ITERATIONS, N, AP, FL, P_LS, false);

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
    cout<<best_contigs<<"\t\t"<<mean_contigs<<"\t\t"<<worst_contigs<<endl;    */
    /////////////////////////////////////////////////////////////////////////////////////////////// 


}