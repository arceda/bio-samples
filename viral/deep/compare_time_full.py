# this compara el tiempo de procesamiento de kameris, castor, mldsp y CNN
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

from numpy.random import seed
seed(1)
import tensorflow as tf
tf.random.set_seed(1)

import pickle
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

from pyblast import BioBlast
from pyblast.utils import make_linear, make_circular
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import glob

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

    ####### para blast ######################################################
    #########################################################################
    queries = [SeqRecord(Seq(sequences[0]))]
    subjects = []
    files_seq = glob.glob( dataset_path + "/" + dataset + "/seq/*" )
    # rega usa 87 instancias, pero en nuestro caso, hay una BD con solo 76 muestras
    for seq_index in range(76):
        #print(files_seq[seq_index])
        seqs_records = SeqIO.parse(files_seq[seq_index], "fasta")

        for record in seqs_records:
            subjects.append( record )
        
    
    t0 = time.time()

    subjects = make_linear(subjects)
    queries = make_linear(queries)

    blast = BioBlast(subjects, queries)
    results = blast.blastn()

    times_blast_acc = time.time() - t0    
    print("times_blast_acc", times_blast_acc)
    #########################################################################
    
    times_castor_acc = 0
    times_kameris_acc = 0
    times_mldsp_acc = 0
    times_cnn1_acc = 0
    times_cnn2_acc = 0
    times_cnn3_acc = 0
    
    # para castor solo
    train = []
    for i in range(sequences.shape[0]):
        train.append([0, sequences[i], 0])   
    k_mers_castor   = fe.generate_K_mers(train, k=5) 
    
    #model_castor = pickle.load(open( 'models/' + dataset + '-castor.joblib', 'rb'))
    #model_kameris = pickle.load(open( 'models/' + dataset + '-kameris.joblib', 'rb'))
    #model_mldsp = pickle.load(open( 'models/' + dataset + '-mldsp.joblib', 'rb'))
    #model_cnn1 = keras.models.load_model('models/' + dataset + '-cnn-tiny.h5', 'rb')
    #model_cnn2 = keras.models.load_model('models/' + dataset + '-cnn-medium.h5', 'rb')
    #model_cnn3 = keras.models.load_model('models/' + dataset + '-cnn-complex.h5', 'rb')

    for i in range(sequences.shape[0]): 
        # castor time        
        t_0 = time.time()
        train_castor = [[0, sequences[i], 0]]   
        X_train, y_train    = fe.generateXYMatrice(train_castor, k_mers_castor, k=5)         
        #pred = model_castor.predict(X_train)
        times_castor_acc += time.time() - t_0

        # kameris time
        t_0 = time.time()
        k_mers_frecuencies = kam.cgr(sequences[i], k=5) 
        #pred = model_kameris.predict( [ k_mers_frecuencies ] )
        times_kameris_acc += time.time() - t_0

        # mldsp time
        t_0 = time.time()
        ns_new, fourier_transform, magnitud_spectra = descriptor(sequences[i], min_seq_len)
        pearsoncorr = np.corrcoef(np.matrix([magnitud_spectra])) # simulacion solo para una instancia, pero se nceesitraria toda la bd
        X = (1 - pearsoncorr)/2
        #pred = model_kameris.predict( [ X ] )
        times_mldsp_acc += time.time() - t_0

        # cnn 1
        t_0 = time.time()
        chaos = chaos_game_representation(probabilities(str(sequences[i]), count_kmers(str(sequences[i]), 5), 5), 5)
        #pred = model_cnn1.predict([chaos])
        times_cnn1_acc += time.time() - t_0

        # cnn 2
        t_0 = time.time()
        chaos = chaos_game_representation(probabilities(str(sequences[i]), count_kmers(str(sequences[i]), 5), 5), 5)
        #pred = model_cnn2.predict([chaos])
        times_cnn2_acc += time.time() - t_0

        # cnn 3
        t_0 = time.time()
        chaos = chaos_game_representation(probabilities(str(sequences[i]), count_kmers(str(sequences[i]), 5), 5), 5)
        #pred = model_cnn3.predict([chaos])
        times_cnn3_acc += time.time() - t_0

    time_kameris = times_kameris_acc/sequences.shape[0]
    time_castor = times_castor_acc/sequences.shape[0]
    time_mldsp = times_mldsp_acc/sequences.shape[0]
    time_cnn1 = times_cnn1_acc/sequences.shape[0]
    time_cnn2 = times_cnn2_acc/sequences.shape[0]
    time_cnn3 = times_cnn3_acc/sequences.shape[0]
    print(times_kameris_acc/sequences.shape[0])
    print(times_castor_acc/sequences.shape[0])
    print(times_mldsp_acc/sequences.shape[0])
    print(times_cnn1_acc/sequences.shape[0])
    print(times_cnn2_acc/sequences.shape[0])
    print(times_cnn3_acc/sequences.shape[0])
    

    with open(current_dir + '/results_v4/time_full.csv', "a") as myfile:    
        myfile.write("\n" + dataset +", " + str(times_blast_acc) + "," + str(time_kameris) + "," + str(time_castor) + ","+ str(time_mldsp) + ","+ str(time_cnn1) + ","+ str(time_cnn2)+ ","+ str(time_cnn3) + "," + str(min_seq_len))

    
