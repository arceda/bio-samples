import numpy as np
import random



def OXFL(crow1, crow2, FL): 
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
 

## PROBADO CON EL EJEMPLO DEL PAPER
if __name__ == "__main__" :
    crow = np.array([9,11,5,8,2,13,14,1,4,3,7,10,6,12])
    victim = np.array([4,10,7,2,13,8,11,5,14,12,6,1,3,9])
    sol = np.zeros(crow.shape[0])
    FL = 0.75

    #C1 y C2 van de 0 a N-1
    C1 = 2
    C2 = 11
    sol = OXFL(crow, victim, FL)
    print("crow", crow)
    print("victim", victim)
    print("sol", sol)

 