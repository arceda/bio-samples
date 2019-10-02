import numpy as np
import random
from local_search import *
from oxfl import OXFL


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
            crows[i] = OXFL(crows[i], crows[random_crow], FL)

            #################     local search     ###################
            r_ls = random.random()
            if r_ls >= P_LS:
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
