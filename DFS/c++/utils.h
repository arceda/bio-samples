#ifndef UTILS_H 
#define UTILS_H

#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <vector>
#include <string>
#include <algorithm>    // copy
#include <iterator>     // ostream_operator

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

#endif