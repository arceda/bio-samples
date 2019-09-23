import numpy as np
import random

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

ITERATIONS = 2
N = 3
AP = 0.02
FL = 0.75
P_LS = 0.49
crows = np.zeros(shape=(N,num_fragments))
memory = np.zeros(shape=(N,num_fragments))

def show_crows(solutions):
    print("SOLUTIONS")
    for i in xrange(solutions.shape[0]):
        print(solutions[i])
        print(fitness(solutions[i]))


def init_population():
    for i in xrange(N):
        crow = np.arange(num_fragments) #individual = [0, 1, 2, ...] each index is a fragment
        np.random.shuffle(crow) #shuffle the fragment, this a ramdon solution
        crows[i] = crow
    return crows


def fitness(solution):
    #print("calculating fitness of: ", solution)
    overlap = 0
    for i in xrange(num_fragments - 1):
        #print("calculating overlap of: ", int(solution[i]), int(solution[i+1]))
        #print("overlap: ", matrix[ int(solution[i]), int(solution[i+1]) ] )
        overlap += matrix[int(solution[i]), int(solution[i+1])] #the overload is yet calculated in matrix
    
    print("fitness calculated: ", overlap)
    return overlap


def OXFL(crow, food):
    n = num_fragments
    C1 = int(random.randint(0, n-1))
    C2 = int((C1+n*FL) % n)
    print("C1, C2: ", C1, C2)

    new_solution = np.zeros(num_fragments)

    if C1+n*FL <= n:
        print("case A OXFL")
        new_solution[C1:C2] = crow[C1:C2]
        #for i in xrange(num_fragments):

        
    else:
        print("case B OXFL")
        new_solution[C1:C2] = food[C1:C2]


    


crows = init_population()
memory = crows.copy()
#print(crows)

i=0
while i < ITERATIONS:
    for i in xrange(N):
        print('CROW ', i)
        random_crow = random.randint(0, N-1) #chose a random crow
        r = random.random()
        if r > AP:
            print("the crow look up", i)
            OXFL(crows[i], random_crow)
        else:
            print("the crow move to ramdon position", i)
            #the crow go to a random position
            #print('the crow go to random position', i, crows[i])
            np.random.shuffle(crows[i])
            #print('the crow went to random position', i, crows[i])
            print("memory[i]: ",i,  memory[i])
            print("crows[i]: ",i,  crows[i])

        if fitness(crows[i]) > fitness(memory[i]):
            print("the new position is better, updating memory")
            memory[i] = crows[i].copy()
            #print("memory[i]: ",i,  memory[i])
            #print("crows[i]: ",i,  crows[i])
        
    i += 1


#show_crows(crows)
#show_crows(memory)
