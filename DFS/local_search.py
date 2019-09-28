import numpy as np

def PALS(N, individual, matrix_w):
    L = []
    for i in range(N-2):
        for j in range(i+1, N-1):
            delta_c, delta_f = calculateDeltas(individual, i, j, matrix_w)
            if delta_c < 0 or (delta_c == 0 and delta_f > 0):
                L.append( (i, j, delta_f, delta_c) )

    if len(L) > 0:
        i, j = selectMovement(L)
        individual = applyMovement(individual, i, j)
    
    return individual

def calculateDeltas(individual, i, j, matrix_w):
    cutoff = 45
    delta_c = delta_f = 0
    delta_f = delta_f - matrix_w[individual[i], individual[i-1]] - matrix_w[individual[j], individual[j+1]]
    delta_f = matrix_w[individual[i-1], individual[j]] - matrix_w[individual[i], individual[j+1]]
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
    print()

def applyMovement(individual, i, j):
    print()