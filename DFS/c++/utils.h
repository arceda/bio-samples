#ifndef UTILS_H 
#define UTILS_H

#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <vector>
#include <string>
#include <algorithm>    // copy
#include <iterator>     // ostream_operator
#include <random>

#include <boost/tokenizer.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;
using namespace boost::numeric; //ublast
using namespace boost::numeric::ublas;//matrix

int CUTOFF = 30;

template<typename T, typename U>
void assign_to_matrix(ublas::matrix<T>& m,std::size_t r,std::size_t c,U const& data)
{
    m.resize(std::max(m.size1(), r+1), std::max(m.size2(), c+1));
    m(r, c) = data;
}

ublas::matrix<int>  read_csv(string file){
    //string data("example.csv");
    ublas::matrix<int> m;

    ifstream in(file.c_str());
    if (!in.is_open()){
        std::cout<<"Could no open csv"<<std::endl;
        return m;
    } 

    typedef tokenizer< escaped_list_separator<char> > Tokenizer;

    std::vector< string > vec;
    string line;
    int i = 0;
    int j = 0;
    while (getline(in,line))
    {
        Tokenizer tok(line);
        vec.assign(tok.begin(),tok.end());

        if (vec.size() < 3) continue;

        //copy(vec.begin(), vec.end(), ostream_iterator<string>(cout, "|"));
        //cout << "\n----------------------" << endl;
        j = 0;
        for (int index = 0; index < vec.size(); index++){
            assign_to_matrix(m, i, j, boost::lexical_cast<int>(vec[index]));             
            j++;
        }        
        i++;
    }

    //std::cout << m << std::endl;
    return m;
}

int consensus(ublas::matrix<int> m, ublas::vector<int> individual){ 
    int contigs = 1;
    for(int i = 0; i < individual.size()-1; i++){
        if ( m(individual[i], individual[i+1]) < CUTOFF )
            contigs++;
    }
    return contigs;
}

int fitness(ublas::matrix<int> m, ublas::vector<int> individual){ 
    int overlap = 0;
    for(int i = 0; i < individual.size()-1; i++){
        overlap += m(individual[i], individual[i+1]);
    }
    return overlap;
}

int min_vector(ublas::vector<int> individual){
    if (individual.size() > 0){
        int min = individual[0];
        int min_index = 0;
        for (int i = 0; i < individual.size(); i++){
            if (individual[i] < min){
                min = individual[i];
                min_index = i;
            }
        }
        return min_index;
    }
    else
        return -1;
}

int max_vector(ublas::vector<int> individual){
    if (individual.size() > 0){
        int max = individual[0];
        int max_index = 0;
        for (int i = 0; i < individual.size(); i++){
            if (individual[i] > max){
                max = individual[i];
                max_index = i;
            }
        }
        return max_index;
    }
    else
        return -1;
}

ublas::vector<int> create_shuffle_vector(int num_fragments){
    std::vector<int> individual(num_fragments);
    for (int i = 0; i < individual.size(); i++)
        individual[i] = i;
    std::random_device rd;
    std::mt19937 randEng(rd());
    std::shuffle(individual.begin(), individual.end(), randEng);

    //from std:vec to ublas::vec
    ublas::vector<int> individual_t(individual.size());
    std::copy(individual.begin(), individual.end(), individual_t.begin());

    return individual_t;
}

int random_number(int low, int high){
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist6(low,high); // distribution in range [1, 6]

    //std::cout << dist6(rng) << std::endl;
    return dist6(rng);
}

float random_float(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    //std::cout << dist6(rng) << std::endl;
    return dis(gen);
}

int find_in_vector(ublas::vector<int> vec, int data){
    for (int i = 0; i < vec.size(); i++){
        if(vec[i] == data)
            return i;
    }
    return -1;
}

void sort_by_column(ublas::matrix<int> &data, int col){

    //cout<<"SORTING ..."<<endl;
    //boost ,atrix to std matrix
    std::vector<std::vector<int> > M;
    for(int i = 0; i < data.size1(); i++){
        std::vector<int> tmp( data.size2());
        for(int j = 0; j < data.size2(); j++){
            tmp[j] = data(i, j);
        }
        M.push_back(tmp);
    }

    /*cout<<"matrix boost"<<data<<endl;
    cout<<"matrix std"<<endl;
    for(int i = 0; i < data.size1(); i++){
        for(int j = 0; j < data.size2(); j++){
            cout<<" "<<M[i][j];
        }
        cout<<endl;
    }*/
   
    std::sort(M.begin(),
              M.end(),
              [col](const std::vector<int>& lhs, const std::vector<int>& rhs) {
                  return lhs[col] > rhs[col];
              });

    
    for(int i = 0; i < data.size1(); i++){
        for(int j = 0; j < data.size2(); j++){
            data(i, j) = M[i][j];
        }
    }
    

}

#endif