import numpy as np

CUTOFF = 30

def consensus(matrix, individual):
        individual = individual.astype(int)
        x = individual.shape[0]

        contigs = 1        
        for i in range(x-1):
                if matrix[individual[i], individual[i+1]] < CUTOFF:
                        contigs += 1

        return contigs

def fitness(matrix, individual):
        individual = individual.astype(int)
        x = individual.shape[0]

        #print("calculating fitness of: ", solution)
        overlap = 0
        for i in range(x-1):
                #print("calculating overlap of: ", int(solution[i]), int(solution[i+1]))
                #print("overlap: ", matrix[ int(solution[i]), int(solution[i+1]) ] )
                overlap += matrix[int(individual[i]), int(individual[i+1])] #the overload is yet calculated in matrix

                #print("fitness calculated: ", overlap)
        return overlap + matrix[int(individual[x-1]), int(individual[0])] 