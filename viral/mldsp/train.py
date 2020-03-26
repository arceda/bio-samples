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
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn import svm
from sklearn.model_selection import train_test_split
import joblib
import os



path_database = '/home/vicente/projects/BIOINFORMATICS/MLDSP/DataBase'
database_name = 'Influenza'
path_database = sys.argv[1]
database_name = sys.argv[2]

def numMappingPP(nucleotide):
    if nucleotide == 'A':
        return -1
    elif nucleotide == 'C':
        return 1
    elif nucleotide == 'G':
        return -1
    elif nucleotide == 'T':
        return 1

# original matlan function
# function [AcNmb, Seq, numberOfClusters, clusterNames, pointsPerCluster] = readFasta(dataSet)
# AcNmb: list of seqId de todas
# Seq: list of sequences in the dataset de todas
# numberOfClusters: numero de clases
# clusterNames: clases
# muestras por clase
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

# Compute the magnitud spectrum of FFT of a numerical sequence
# seq: sequence, string of letters A, C, G, T
# median_len: median size of all sequences in a dataset
def descriptor(seq, median_len):
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

    return  ns_new, fourier_transform, magnitud_spectra


if __name__ == "__main__" :
    sequences, number_of_clases, cluster_names, points_per_cluster = readFasta(database_name)

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
        ns_new, fourier_transform, magnitud_spectra = descriptor(seq, median_len)
        
        nm_val_SH.append(ns_new)
        f.append(fft(fourier_transform))
        lg.append(magnitud_spectra)
        

    #################################################################################################
    # distance calculation by Pearson correlation coefficient
    print('Computing Distance matrix .... ...')

    # pandas version
    #lg_df = pd.DataFrame(np.transpose(lg)) # transpose in order to compute PCC by observation
    #pearsoncorr = lg_df.corr(method='pearson')  #Pearson correlation coefficient [-1 1]
    #dist_mat = (1 - pearsoncorr)/2  #  normalize between 0 and 1

    # numpy version
    pearsoncorr = np.corrcoef(np.matrix(lg))
    dist_mat = (1 - pearsoncorr)/2

    # testing Pearson Correlation Distance
    #X = [[1460, 517.201, 230.163, 453.649, 266.169, 267.257], [1340, 569.351, 219.907, 473.615, 239.587, 229.557], [1462, 622.617, 324.276, 503.927, 214.432,223.652], [1994, 672.012, 456.685, 412.211, 219.068, 131.52]]
    #X_pd = pd.DataFrame(np.transpose(np.matrix(X)))
    #X_corr = X_pd.corr(method='pearson')
    #X_corr_norm = (1 - X_corr)/2

    #X_corr_2 = np.corrcoef(np.matrix(X))
    #X_corr_norm_2 = (1 - X_corr_2)/2
    ##################################################################################################
    # train

    kf = KFold(n_splits=5)
    X = dist_mat
    y = sequences[:, 2]

    clf = svm.SVC(kernel='linear', C=1)
    scores = cross_val_score(clf, X, y, cv=5)
    print("scores cv=5", scores)
    print("mean score", statistics.mean(scores))

    # thain the whole database
    clf = svm.SVC(kernel='linear', C=1)
    clf.fit(X, y)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    file_name = current_dir + '/models/' + database_name
    joblib.dump(clf, file_name + ".sav")
    np.savetxt(file_name + "_magnitud_spectrum.csv", lg, delimiter=',', fmt='%f')

    print(len(lg), " features")

    '''
    for train_index, test_index in kf.split(dist_mat):
        print('TRAIN:', train_index, 'TEST:', test_index)
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        
    '''
    