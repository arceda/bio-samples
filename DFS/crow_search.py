import numpy as np
import random


#range es mas eficiente que range in python2, in python3 range = range()pyton2

###########################################################################################
############################### read data ##############################################
instance = 'x60189_7'



matrix = np.genfromtxt(instance + '/matrix_conservative.csv', delimiter=',')
fragments = np.array([])

fh = open(instance + '/frag_x60189_7.dat')
while True:
    line = fh.readline().replace('\n', '')
    fragments = np.append(fragments, line)
    if not line:
        break
fh.close()
fragments = fragments[1::2]
num_fragments = fragments.shape[0]

#print(matrix.shape)
#print(fragments.shape)

###########################################################################################
###########################################################################################

ITERATIONS = 300
N = 32
AP = 0.02
FL = 0.75
P_LS = 0.49
crows = np.zeros(shape=(N,num_fragments))
memory = np.zeros(shape=(N,num_fragments))

def show_crows(solutions):
    print("SOLUTIONS")
    for i in range(solutions.shape[0]):
        print(solutions[i])
        print(fitness(solutions[i]))


def init_population():
    for i in range(N):
        crow = np.arange(num_fragments) #individual = [0, 1, 2, ...] each index is a fragment
        np.random.shuffle(crow) #shuffle the fragment, this a ramdon solution
        crows[i] = crow
    return crows


def fitness(solution):
    #print("calculating fitness of: ", solution)
    overlap = 0
    for i in range(num_fragments - 1):
        #print("calculating overlap of: ", int(solution[i]), int(solution[i+1]))
        #print("overlap: ", matrix[ int(solution[i]), int(solution[i+1]) ] )
        overlap += matrix[int(solution[i]), int(solution[i+1])] #the overload is yet calculated in matrix
    
    #print("fitness calculated: ", overlap)
    return overlap


def OXFL(crow1, crow2):
    crow = crow1.copy()
    victim = crow2.copy()

    n = num_fragments
    #C1 y C2 van de 0 a num_fragments-1
    C1 = int(random.randint(0, n-1))
    C2 = int((C1+(n-1)*FL) % (n-1))
    
    #print("C1, C2: ", C1, C2)

    sol = np.zeros(num_fragments)

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

def P2M_F(individual):
    print("local search")


crows = init_population()
memory = crows.copy()
#print(crows)

iter=0
while iter < ITERATIONS:
    print("ITERATION: ", iter)
    for i in range(N):
        #print('CROW ', i)
        random_crow = random.randint(0, N-1) #chose a random crow
        r = random.random()
        if r >= AP:
            #print("the crow look up", i)
            #print("perform oxfl operator")
            crows[i] = OXFL(crows[i], crows[random_crow])

            #################     local search     ###################
            r_ls = random.random()
            if r_ls >= P_LS:
                P2M_F(crows[i])

        else:
            #print("the crow move to ramdon position", i)
            #the crow go to a random position
            #print('the crow go to random position', i, crows[i])
            np.random.shuffle(crows[i])
            #print('the crow went to random position', i, crows[i])
            #print("memory[i]: ",i,  memory[i])
            #print("crows[i]: ",i,  crows[i])

        if fitness(crows[i]) > fitness(memory[i]):
            #print("the new position is better, updating memory")
            memory[i] = crows[i].copy()
            #print("memory[i]: ",i,  memory[i])
            #print("crows[i]: ",i,  crows[i])
        
    iter += 1


# get the best fitness
best_fitness = 0
for i in range(N):
    fit = fitness(memory[i])
    if fit > best_fitness:
        best_fitness = fit

print("best fitness", best_fitness)
