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

from train import descriptor

database_name = "Primates"
path_seq = sys.argv[1]
database_name = sys.argv[2]


current_dir = os.path.dirname(os.path.abspath(__file__))
file_name = current_dir + '/models/' + database_name

print("readding magnitud spectrum...")
lg = np.genfromtxt(file_name + "_magnitud_spectrum.csv", delimiter=',', dtype=float)
print(lg.shape)
median_len = lg.shape[1]

sequences = []
seqs = SeqIO.parse(path_seq, "fasta")                       
for record in seqs:
    sequences.append([record.id, record.seq.upper()])     

seq_distances = []              
for seq in sequences:
    #print(seq[1])
    ns_new, fourier_transform, magnitud_spectra = descriptor(seq[1], median_len)
    #print(magnitud_spectra, len(magnitud_spectra))
    distances = [] # the feature vector
    for observation in lg:
        #print(observation.shape)
        r = np.corrcoef(observation, magnitud_spectra)[0, 1] #corrcoef return a matrix
        distances.append((1-r)/2)

    seq_distances.append(distances)

print(clf.predict(seq_distances))
    

