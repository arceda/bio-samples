import numpy as np

crow = np.array([9,11,5,8,2,13,14,1,4,3,7,10,6,12])
victim = np.array([4,10,7,2,13,8,11,5,14,12,6,1,3,9])
sol = np.zeros(14)

N = 3
AP = 0.02
FL = 0.75
P_LS = 0.49
num_fragments = 14

def OXFL(crow, victim):
    n = num_fragments
    #C1 y C2 van de 0 a num_fragments-1
    C1 = int(random.randint(0, n-1))
    C2 = int((C1+n*FL) % n)
    print("C1, C2: ", C1, C2)

    new_solution = np.zeros(num_fragments)

    if C1+n*FL <= n:
        print("case A OXFL")
        X = crow
        M = victim   
    else:
        print("case B OXFL")
        X = victim
        M = crow

    # copy from C1 to C2 to new solution
    sol[C1:C2] = X[C1:C2]
    
    #delete the range C1:c2 from M
    for i, x in enumerate(X[C1:C2]):
        index, = np.where(M == x)
        M = np.delete(M, index)  

    #print("M", M)
    sol[0:C1] = M[0:C1] # insert the firsts to C1 to new solution
    M = np.delete(M, range(C1)) # delete the elements inserted
    #print("M", M)
    sol[C2:sol.shape[0]] = M
    #print("sol", sol)
 

## PROBADO CON EL EJEMPLO DEL PAPER
#C1 y C2 van de 0 a N-1
C1 = 2
C2 = 11

if True:
    print("case A OXFL")
    X = crow
    M = victim      
else:
    print("case B OXFL")
    X = victim
    M = crow

# copy from C1 to C2 to new solution
sol[C1:C2] = X[C1:C2]

print("X", X)
print("M", M)
print("sol", sol)

#delete the range C1:c2 from M
for i, x in enumerate(X[C1:C2]):
    index, = np.where(M == x)
    M = np.delete(M, index)
    #print(index, M[index])
 

print("M", M)
sol[0:C1] = M[0:C1] # insert the firsts to C1 to new solution
M = np.delete(M, range(C1)) # delete the elements inserted
print("M", M)
sol[C2:sol.shape[0]] = M

print("sol", sol)
 