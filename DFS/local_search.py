import numpy as np
from utils import fitness
from utils import consensus

# IMPLEMENTACION DE PALS CON SUS MODIFICACIONES
######################################### REFERENCES ##############################################
#[1] "A New Local Search Algorithm for the DNA Fragment Assembly Problem"
#[2] "An improved problem aware local search algorithm for the DNA fragment assembly problem"
#[3] "A hybrid crow search algorithm for solving the DNA fragment assembly problem"
###################################################################################################


def PALS(K, individual, matrix_w):
    #print("individual receive in PALS: ", individual)
    individual = individual.astype(int)

    iterations = 0
    #while iterations < 300:
    
    while iterations < 3000:
        
        L = []
        for i in range(1, K):
            for j in range(0, K-1):
                delta_c, delta_f = calculateDeltas(individual, i, j, matrix_w)
                
                ###################################################################################################
                #if delta_c < 0 or (delta_c == 0 and delta_f > 0): #PALS original [1, 2]
                if delta_f > 0: #PALS modificado en [3]
                    L.append( [i, j, delta_f, delta_c] )

        if len(L) > 0:
            ###################################################################################################
            #PALS original [1]
            #i, j = selectMovement(L)
            #individual = applyMovement(individual, i, j)

            ###################################################################################################
            #PALS modificado [2]
            #individual = applyMovement_PALS2many(individual, L)

            ###################################################################################################
            #PALS modificado en [3]
            individual = applyMovement_PALS2many_fit(individual, L)
            #break

        #print(" interation PALS: ", iterations, " candidates number: ", len(L), " fitness: ", fitness(matrix_w, individual))

        iterations += 1

        if len(L) <= 0:
            break
    
    return individual

def calculateDeltas(individual, i, j, matrix_w):
    cutoff = 30
    delta_c = delta_f = 0
    delta_f = delta_f - matrix_w[individual[i-1], individual[i]] - matrix_w[individual[j], individual[j+1]]
    delta_f = delta_f + matrix_w[individual[i-1], individual[j]] + matrix_w[individual[i], individual[j+1]]
    if matrix_w[individual[i-1], individual[i]] > cutoff:
        delta_c = delta_c + 1
    if matrix_w[individual[j], individual[j+1]] > cutoff:
        delta_c = delta_c + 1
    if matrix_w[individual[i-1], individual[j]] > cutoff:
        delta_c = delta_c - 1
    if matrix_w[individual[i], individual[j+1]] > cutoff:
        delta_c = delta_c - 1
    
    return delta_c, delta_f

def selectMovement(L):
    # get the posible movement with minimun delta_c
    x = len(L)
    L_temp = np.matrix(L)
    delta_c_list = L_temp[:,3]
    min_delta_c = np.amin(delta_c_list)

    L_with_min_delta_c = []
    for i in range(x):
        if L_temp[i,3] == min_delta_c:
            L_with_min_delta_c.append(np.squeeze(np.asarray(L_temp[i,:])))

       
    # get the posible movement with maximun delta_f
    x = len(L_with_min_delta_c)
    L_temp = np.matrix(L_with_min_delta_c)
    delta_f_list = L_temp[:,2]
    max_delta_f = np.amin(delta_f_list)
    
    L_with_max_delta_f = []
    for i in range(x):
        if L_temp[i,2] == max_delta_f:
            L_with_max_delta_f.append(np.squeeze(np.asarray(L_temp[i,:])))

    L_temp = np.matrix(L_with_max_delta_f)
    
    #print(L_temp.shape)   
    #print(L_temp)
    return int(L_temp[0, 0]), int(L_temp[0, 1])
    

def applyMovement(individual, i, j):
    #print("applying movement")
    #print(i, j)
    #print(individual)
    tmp = individual[i]
    individual[i] = individual[j]
    individual[j] = tmp
    #print(individual)
    return individual

# aplica el algoritmo de PALS2many* propuesto en [2]
def applyMovement_PALS2many(individual, L):
    L = np.array(L)
    L = L.astype(int)    
    #np.savetxt("L_initial.csv", L, delimiter=",", fmt='%i')

    #sort ascending by delta_c (4th column), them by delta_f (3th column)
    L = L[np.lexsort((L[:,2], L[:,3]))]    
    #np.savetxt("L_final.csv", L, delimiter=",", fmt='%i')

    #applied all the movement according to "An improved problem aware local search algorithm for the DNA fragment assembly problem"
    index_used = np.array([]) # save the movements applied
    count=0
    for posible_movement in range(L.shape[0]):
        #print("movement:", L[posible_movement])
        #print("individual befores movement:", individual)
        i = L[posible_movement][0]
        j = L[posible_movement][1]

        #si el movieminto no fue aplicado antes
        #print(np.where(index_used==i)[0], np.where(index_used==j)[0])
        if len(np.where(index_used==i)[0]) == 0 and len(np.where(index_used==j)[0]) == 0:
            #print("apply movement")
            index_used = np.append(index_used, i)
            index_used = np.append(index_used, j)
            
            #swap
            tmp = individual[i]
            individual[i] = individual[j]
            individual[j] = tmp
            count += 1
            #print("individual after movement:", individual)

   
    #print("movements applied: ", count)
    #print("new individual: ", individual) 
    return individual

# aplica el algoritmo de PALS2many*fit propuesto [3]
def applyMovement_PALS2many_fit(individual, L):
    L = np.array(L)
    L = L.astype(int)    
    #np.savetxt("L_initial.csv", L, delimiter=",", fmt='%i')

    #sorting descending by delta_f
    L = L[L[:,2].argsort()[::-1]]   
    #np.savetxt("L_final.csv", L, delimiter=",", fmt='%i')

    #applied all the movement according to "An improved problem aware local search algorithm for the DNA fragment assembly problem"
    index_used = np.array([]) # save the movements applied
    count=0
    for posible_movement in range(L.shape[0]):
        #print("movement:", L[posible_movement])
        #print("individual befores movement:", individual)
        i = L[posible_movement][0]
        j = L[posible_movement][1]

        #si el movieminto no fue aplicado antes
        #print(np.where(index_used==i)[0], np.where(index_used==j)[0])
        if len(np.where(index_used==i)[0]) == 0 and len(np.where(index_used==j)[0]) == 0:
            #print("apply movement")
            index_used = np.append(index_used, i)
            index_used = np.append(index_used, j)
            
            #swap
            tmp = individual[i]
            individual[i] = individual[j]
            individual[j] = tmp
            count += 1
            #print("individual after movement:", individual)

   
    #print("movements applied: ", count)
    #print("new individual: ", individual) 
    return individual

if __name__ == "__main__" :
    instance = 'x60189_4'
    matrix = np.genfromtxt(instance + '/matrix_conservative.csv', delimiter=',')
    #print(matrix.shape)
    num_fragments = matrix.shape[0]
    aleatory_solution = np.arange(num_fragments)
    np.random.shuffle(aleatory_solution)

    print("initial solution: ", aleatory_solution)    
    sol = PALS(num_fragments, aleatory_solution, matrix)
    fitness_temp = fitness(matrix, sol)
    contigs_temp = consensus(matrix, sol)
    print("fitness: ", fitness_temp, "contig: ", contigs_temp)
    print("final solution: ", sol)

    #################################################### pruebas #########################################
    num_test = 30.0
    fitness_acum = best_fitness = contig_acum = 0.0
    best_contig = 100

    """
    print("TESTING 30 ITEARTIONS...")    

    for i in range(int(num_test)):
        aleatory_solution = np.arange(num_fragments)
        np.random.shuffle(aleatory_solution)

        print("testing...", i)
        sol = PALS(num_fragments, aleatory_solution, matrix)
        fitness_temp = fitness(matrix, sol)
        contigs_temp = consensus(matrix, sol)

        print("fitness: ", fitness_temp, "contig: ", contigs_temp)

        fitness_acum += fitness_temp
        contig_acum += contigs_temp
        if fitness_temp > best_fitness:
            best_fitness = fitness_temp
        if contigs_temp < best_contig:
            best_contig = contigs_temp

    fitness_mean = fitness_acum/num_test
    contigs_mean = contig_acum/num_test    
    
    print("fitness_mean: ", fitness_mean)
    print("contigs_mean: ", contigs_mean)"""
    