# MLDSP: machine learning with DIgiral Signal Processing
# A python implementation of MLDSP

import numpy as np
from Bio import SeqIO
import re
import sys
import glob
import statistics
import pywt
#from sympy import fft es muy lento
from scipy.fftpack import fft
import pandas as pd

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

    sequences_mat = np.array(sequences)

    return sequences_mat, number_of_clases, cluster_names, points_per_cluster


sequences, number_of_clases, cluster_names, points_per_cluster = readFasta('Primates')

#calculate length stats
sequences_size  = list(map(len, sequences[:, 1]))
total_seq       = len(sequences_size)

max_len     = max(sequences_size)
min_len     = min(sequences_size)
mean_len    = int(statistics.mean(sequences_size))
median_len  = int(statistics.median(sequences_size))
#print("max_len, min_len, mean_len, median_len ", max_len, min_len, mean_len, median_len)

nm_val_SH     = []
f           = []
lg          = []

print('Generating numerical sequences, applying DFT, computing magnitude spectra ...')

for seq in sequences[:, 1]:

    ns  = list(map(numMappingPP, seq))  # here we map the nucleotides to numbers

    ## we ensure that all the sequences have the same size 
    I   = median_len - len(ns) # change "median_len" to other length stat for length normalization
    if I > 0: # when the seq size is smaller than the median
        ns_temp = pywt.pad(ns, I, 'antisymmetric') #wextend('1','asym',ns,I);
        ns_temp = np.array(ns_temp)
        #print(len(ns_temp), len(ns))
        ns_new = ns_temp[I:ns_temp.shape[0]] #nsNew = nsTemp((I+1):length(nsTemp));        
    elif I < 0: # when the seq size is bigger than the median
        ns = np.array(ns)
        ns_new = ns[0:median_len]
    else:
        ns_new = np.array(ns)

    #print(ns_new.shape)
    #print("processing FFT...")
    fourier_transform = fft(ns_new)
    #print("getting  magnitude spectra...")
    magnitud_spectra = np.abs(fourier_transform) # %magnitude spectra

    #print(len(fourier_transform))
    #print(len(magnitud_spectra))

    nm_val_SH.append(ns_new)
    f.append(fft(fourier_transform))
    lg.append(magnitud_spectra)
    

#################################################################################################
# distance calculation by Pearson correlation coefficient
print('Computing Distance matrix .... ...')

lg_df = pd.DataFrame(lg)
pearsoncorr = lg_df.corr(method='pearson')
print(pearsoncorr.shape)