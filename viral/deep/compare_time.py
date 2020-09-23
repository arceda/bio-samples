# this compara el tiempo de procesamiento de kameris, castor, mldsp y CNN
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

from numpy.random import seed
seed(1)
import tensorflow as tf
tf.random.set_seed(1)


from sklearn.model_selection import KFold 
from sklearn.model_selection import train_test_split
import collections
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm
import pylab
import math
import os
import sys
import cv2
import numpy as np

from Bio import SeqIO
import re
from cgr import save_fcgr
import glob


from tensorflow.keras import datasets, layers, models
from keras.models import Sequential
from keras.layers import Conv2D, MaxPooling2D, Dense, Flatten, Dropout
from keras.utils import to_categorical
from keras.utils.vis_utils import plot_model
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import accuracy_score

from keras.utils import to_categorical
from keras.layers.normalization import BatchNormalization
import seaborn as sns
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import confusion_matrix
from numpy import argmax
import matplotlib.pyplot as plt

import random
import statistics
import pywt
import csv

import mykameris as kam
import feature_extractor as fe
import statistics
import pywt
from scipy.fftpack import fft

from cgr import chaos_game_representation
from cgr import probabilities
from cgr import count_kmers

import time

def numMappingPP(nucleotide):
    if nucleotide == 'A':
        return -1
    elif nucleotide == 'C':
        return 1
    elif nucleotide == 'G':
        return -1
    elif nucleotide == 'T':
        return 1


# Compute the magnitud spectrum of FFT of a numerical sequence
# seq: sequence, string of letters A, C, G, T
# median_len: median size of all sequences in a dataset
def descriptor(seq, median_len):
    #print(seq, median_len) 
    ns  = list(map(numMappingPP, seq))  # here we map the nucleotides to numbers
    ns = list(filter(None.__ne__, ns))

    
    ns_new = np.array(ns)

    #print(ns_new.shape)
    #print("processing FFT...")
    fourier_transform = fft(ns_new)
    #print("getting  magnitude spectra...")
    magnitud_spectra = np.abs(fourier_transform) # %magnitude spectra

    return  ns_new, fourier_transform, magnitud_spectra

def read_seq(path_database, database):
    path = path_database + '/' + database
    X = []
    X_train = []
    y_train = []
    X_test = []
    y_test = []
    #119 train, 29 test

    labels_train = np.loadtxt(path + '/train_labels.csv', dtype=(str,str), delimiter=',')
    labels_test = np.loadtxt(path + '/test_labels.csv', dtype=(str,str), delimiter=',')
    file_labels = np.vstack( (labels_train, labels_test) )
    labels = np.hstack( (labels_train[:,1], labels_test[:,1]) ) #hstack, xq tiene una dimension

    #integer_encoded, onehot_encoded, label_encoder = one_hot_encode(labels)
    
    for file_name, label in file_labels:
        file_name_without_ext = file_name.split('.')[0]
        seqs = SeqIO.parse(path + '/seq/' + file_name, "fasta") 
        
        for record in seqs:  
            X.append( [ str(record.seq.upper()), len( str(record.seq.upper())  ) ] )
            break #read one sequence by file.. all the files have one seq
    
        
    # solo consideramos 10 sequencias
    X = np.array(X)
    np.random.shuffle(X)
    X = X[0:10, :] # solo tomamos 10 secuencias
    X = X[X[:,1].argsort()] # ordenamos por longitud
    min_seq_len = X[0,1]

    # hacemos que todas las secuencias tengan la misma longitud
    for i in range( X.shape[0]):
        X[i, 0] = (X[i, 0])[0:int(min_seq_len)]

    return X[:, 0], min_seq_len


    
    
    #return X_train, y_train, X_test, y_test, set(labels)


dataset_path = sys.argv[1]
#python3 compare_time.py '/home/vicente/datasets/MLDSP/' 

datasets = [    'Primates', 'Dengue', 'Protists', 'Fungi', 'Plants', 'Amphibians', 'Insects', '3classes', 'Vertebrates',
                'HIVGRPCG', 'HIVSUBCG', 'HIVSUBPOL', 'INFSUBHA', 'INFSUBMP', 'INSUBFNA', 'EBOSPECG', 'HBVGENCG', 'RHISPECG', 'HPVGENCG']


current_dir = os.path.dirname(os.path.abspath(__file__))

for i, dataset in enumerate(datasets):

    print(i, "\n\nEVALUATING DATASET: ", dataset, "...")
    print('===========================================================')
    print('===========================================================') 

    sequences, min_seq_len = read_seq(dataset_path, dataset)

    times_castor_acc = 0
    times_kameris_acc = 0
    times_mldsp_acc = 0
    times_cnn_acc = 0
    # para castror solo
    train = []
    for i in range(sequences.shape[0]):
        train.append([0, sequences[i], 0])   
    k_mers_castor   = fe.generate_K_mers(train, k=5) 

    # castor time
    t_0 = time.time()
    X_train, y_train    = fe.generateXYMatrice(train, k_mers_castor, k=5) 
    times_castor_acc += time.time() - t_0

    for i in range(sequences.shape[0]):
        # kameris time
        t_0 = time.time()
        k_mers_frecuencies = kam.cgr(sequences[i], k=5) 
        times_kameris_acc += time.time() - t_0

        # mldsp time
        t_0 = time.time()
        ns_new, fourier_transform, magnitud_spectra = descriptor(sequences[i], min_seq_len)
        times_mldsp_acc += time.time() - t_0

        # cnn
        t_0 = time.time()
        chaos = chaos_game_representation(probabilities(str(sequences[i]), count_kmers(str(sequences[i]), 5), 5), 5)
        times_cnn_acc += time.time() - t_0

    time_kameris = times_kameris_acc/sequences.shape[0]
    time_castor = times_castor_acc/sequences.shape[0]
    time_mldsp = times_mldsp_acc/sequences.shape[0]
    time_cnn = times_cnn_acc/sequences.shape[0]
    print(times_kameris_acc/sequences.shape[0])
    print(times_castor_acc/sequences.shape[0])
    print(times_mldsp_acc/sequences.shape[0])
    print(times_cnn_acc/sequences.shape[0])
    

    with open(current_dir + '/results_v3/time.csv', "a") as myfile:    
        myfile.write("\n" + dataset +", " + str(time_kameris) + "," + str(time_castor) + ","+ str(time_mldsp) + ","+ str(time_cnn) + "," + str(min_seq_len))