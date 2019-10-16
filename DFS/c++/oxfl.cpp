#include <iostream>     // cout, endl
#include <vector>
#include <string>
#include <algorithm>    // copy
#include <iterator>     // ostream_operator
#include <random>

#include <boost/tokenizer.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm/find.hpp>
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

ublas::vector<int> operator_oxfl(ublas::vector<int> crow, ublas::vector<int> victim, float FL){
    
    int n = crow.size();

    //C1 y C2 van de 0 a num_fragments-1
    //int C1 = int(random_number(0, n-1))
    //int C2 = int((C1+(n-1)*FL) % (n-1))
    int C1 = 2;
    int C2 = 11; 
    int temp;

    ublas::vector<int> sol(n);
    ublas::vector<int> X;
    ublas::vector<int> M;

    if (C1 < C2){
        X = crow;
        M = victim;   
    }else{
        X = victim;
        M = crow;
        temp = C1;
        C1 = C2;
        C2 = temp;
    }

    //sol[C1:C2] = X[C1:C2]
    for (int i = C1; i < C2; i++){
        sol[i] = X[i];
    }

    //delete the range C1:c2 from M
    for (int i = C1; i < C2; i++){
        int index = find_in_vector(M, X[i]);
        M.erase_element(index);
        M.resize(M.size()-1);
    }

    /*
    for i, x in enumerate(X[C1:C2]):
        index, = np.where(M == x)
        M = np.delete(M, index)  
*/
    cout<<"sol:"<<sol<<endl;
    cout<<"M:"<<M<<endl;

    return sol;

    /*
    #print("num_fragmetns calculated in oxlf",crow1.shape[0])   
    crow = crow1.copy()
    victim = crow2.copy()

    n = crow1.shape[0]
    #C1 y C2 van de 0 a num_fragments-1
    C1 = int(random.randint(0, n-1))
    C2 = int((C1+(n-1)*FL) % (n-1))
    
    #print("C1, C2: ", C1, C2)

    sol = np.zeros(n)

    #if C1+n*FL <= n-1:
    if C1 < C2:
        #print("case A OXFL")
        X = crow
        M = victim   
    else:
        #print("case B OXFL")
        X = victim
        M = crow
        temp = C1
        C1 = C2
        C2 = temp

    # copy from C1 to C2 to new solution
    #print("X[C1:C2]", X[C1:C2].shape, X[C1:C2])

    sol[C1:C2] = X[C1:C2]
    
    #print("X", X.shape, X)
    #print("M", M.shape, M)
    #print("sol", sol.shape, sol)      
    #delete the range C1:c2 from M
    for i, x in enumerate(X[C1:C2]):
        index, = np.where(M == x)
        M = np.delete(M, index)  

    #print("M", M.shape, M)
    sol[0:C1] = M[0:C1] # insert the firsts to C1 to new solution
    M = np.delete(M, range(C1)) # delete the elements inserted
    #print("M", M.shape, M)
    sol[C2:sol.shape[0]] = M
    #print("sol", sol.shape, sol)    

    return sol
    */

}

int main(){

    std::vector<int> tmp{9,11,5,8,2,13,14,1,4,3,7,10,6,12};
    ublas::vector<int> crow(tmp.size());
    std::copy(tmp.begin(), tmp.end(), crow.begin());

    std::vector<int> tmp2{4,10,7,2,13,8,11,5,14,12,6,1,3,9};
    ublas::vector<int> victim(tmp2.size());
    std::copy(tmp2.begin(), tmp2.end(), victim.begin());
   
    float FL = 0.75;
    
    ublas::vector<int> sol = operator_oxfl(crow, victim, FL);    
    
  
}