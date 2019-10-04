import numpy as np

def PALS(K, individual, matrix_w):
    #print("individual receive in PALS: ", individual)
    individual = individual.astype(int)

    iterations = 0
    while iterations < 300:
        L = []
        for i in range(K-2):
            for j in range(i+1, K-1):
                delta_c, delta_f = calculateDeltas(individual, i, j, matrix_w)
                if delta_c < 0 or (delta_c == 0 and delta_f > 0):
                    L.append( (i, j, delta_f, delta_c) )

        if len(L) > 0:
            i, j = selectMovement(L)
            individual = applyMovement(individual, i, j)

        iterations += 1
    
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
    max_delta_f = np.amax(delta_f_list)
    
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


def fitness(solution, matrix):
    #print("calculating fitness of: ", solution)
    overlap = 0
    for i in range(num_fragments - 1):
        #print("calculating overlap of: ", int(solution[i]), int(solution[i+1]))
        #print("overlap: ", matrix[ int(solution[i]), int(solution[i+1]) ] )
        overlap += matrix[int(solution[i]), int(solution[i+1])] #the overload is yet calculated in matrix
    
    #print("fitness calculated: ", overlap)
    return overlap


if __name__ == "__main__" :
    instance = 'x60189_7'
    matrix = np.genfromtxt(instance + '/matrix_conservative.csv', delimiter=',')
    num_fragments = 68
    aleatory_solution = np.arange(num_fragments)
    np.random.shuffle(aleatory_solution)
    sol = PALS(num_fragments, aleatory_solution, matrix)
    print(fitness(sol, matrix))