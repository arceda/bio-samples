# MLDSP: machine learning with DIgiral Signal Processing
# A python implementation of MLDSP

import numpy as np
from Bio import SeqIO
import re
import sys
import glob
import statistics

# original matlan function
# function [AcNmb, Seq, numberOfClusters, clusterNames, pointsPerCluster] = readFasta(dataSet)
# AcNmb: list of seqId de todas
# Seq: list of sequences in the dataset de todas
# numberOfClusters: numero de clases
# clusterNames: clases
# muestras por clase

path_database = '/home/vicente/projects/BIOINFORMATICS/MLDSP/DataBase'
#path_database = sys.argv[1]

def numMappingPP(nucleotide):
    if nucleotide == 'A':
        return -1
    elif nucleotide == 'C':
        return 1
    elif nucleotide == 'G':
        return -1
    elif nucleotide == 'T':
        return 1


def readFasta(database):
    path = path_database + '/' + database

    number_of_clases = 0
    cluster_names = []
    points_per_cluster = []
    sequences = []
    
    #print(glob.glob(path + '/*' ))
    clusters = glob.glob(path + '/*' )
    number_of_clases = len(clusters)
    for cluster in clusters:
        cluster_name = cluster.split('/')[-1]
        cluster_names.append( cluster_name )

        # read each fasta file
        files = clusters = glob.glob(cluster + '/*.txt' )
        points_per_cluster.append(len(files))
        
        # read sequences from each file, the majority have one sequence per file
          
        for file in files:        
            seqs = SeqIO.parse(file, "fasta") 
                      
            for record in seqs:
                sequences.append([record.id, record.seq.upper(), cluster_name])  
                  
    
    #print(number_of_clases)
    #print(cluster_names)
    #print(points_per_cluster)
    #print(len(sequences))
    #print(sequences)
    sequences_mat = np.array(sequences)

    return sequences_mat, number_of_clases, cluster_names, points_per_cluster

sequences, number_of_clases, cluster_names, points_per_cluster = readFasta('Primates')

#calculate length stats
sequences_size  = list(map(len, sequences[:, 1]))
total_seq       = len(sequences_size)

max_len     = max(sequences_size)
min_len     = min(sequences_size)
mean_len    = statistics.mean(sequences_size)
median_len  = statistics.median(sequences_size)
#print(max_len, min_len, mean_len, median_len)

nmValSH     = np.zeros((total_seq))
f           = np.zeros((total_seq))
lg          = np.zeros((total_seq))

for seq in sequences[:, 1]:
    #print(seq)
    ns  = list(map(numMappingPP, seq))  # here we map the nucleotides to numbers
    I   = median_len - len(ns) # change "median_len" to other length stat for length normalization
    if I > 0:


