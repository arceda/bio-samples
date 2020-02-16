# this script will evaluate the k-mers, number of features an fi score for 100 iteration cooresponding to fig 1, 2 1n 3 
# of paper: Toward an Alignment-Free Method for Feature Extraction and Accurate Classification of Viral Sequences


import csv
import numpy as np
from Bio import SeqIO
import re
import sys

import matplotlib.pyplot as plt
from sklearn.metrics import f1_score
from sklearn.feature_selection import RFE
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold

from sklearn.svm import SVC
from sklearn.naive_bayes import MultinomialNB
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier

import time
import feature_extractor as fe

# at first we will evaluate just the HIV datasets
datasets = ['HIV/HIVGRPCG', 'HIV/HIVSUBCG', 'HIV/HIVSUBPOL'] 
#datasets = ['PAPILLOMA/HPVGENCG', 'PAPILLOMA/HPVSPECG', 'HIV/HIVGRPCG', 'HIV/HIVSUBCG', 'HIV/HIVSUBPOL', 'HEPATITIS-B/HBVGENCG''POLYOMAVIRUS/POLSPEVP1', 'POLYOMAVIRUS/POLSPEVP2', 'POLYOMAVIRUS/POLSPEVP3', 'POLYOMAVIRUS/POLSPEST', 'POLYOMAVIRUS/POLSPELT', 'RHINOVIRUS/RHISPECG', 'DENGE/DENSPECG', 'INFLUENZA/INSUBFNA', 'INFLUENZA/INFSUBHA', 'INFLUENZA/INFSUBMP', 'EBOLA/EBOSPECG']

iterations_number = 2
#dataset_path = "/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/"
dataset_path = sys.argv[1]
f = open("results_k-mers_features.txt", "a")
f.write("Dataset,k,nfeatures,k-mers")


for i, dataset in enumerate(datasets):
    print(i, "EVALUATING DATASET: ", dataset, "...")
    
    for j in range(iterations_number):
        print("iteration ...", j)
        best_k_mers, best_k_length = fe.getBestKmersAndFeatures(dataset_path + dataset)
        print(dataset + "," + str(best_k_length) + "," + str(len(best_k_mers)) + "," + str(best_k_mers))
        f.write(dataset + "," + str(best_k_length) + "," + str(len(best_k_mers)) + "," + str(best_k_mers))

f.close()