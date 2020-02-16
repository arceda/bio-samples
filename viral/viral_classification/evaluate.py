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
from random import shuffle


if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

# at first we will evaluate just the HIV datasets
#datasets = ['POLYOMAVIRUS/POLSPEVP1', 'POLYOMAVIRUS/POLSPEVP2', 'POLYOMAVIRUS/POLSPEVP3'] 
datasets = ['POLYOMAVIRUS/POLSPEST', 'POLYOMAVIRUS/POLSPELT', 'HEPATITIS-B/HBVGENCG', 'HIV/HIVGRPCG', 'HIV/HIVSUBCG', 'HIV/HIVSUBPOL'] 
#datasets = ['PAPILLOMA/HPVGENCG', 'PAPILLOMA/HPVSPECG', 'HIV/HIVGRPCG', 'HIV/HIVSUBCG', 'HIV/HIVSUBPOL', 'HEPATITIS-B/HBVGENCG', 'POLYOMAVIRUS/POLSPEVP1', 'POLYOMAVIRUS/POLSPEVP2', 'POLYOMAVIRUS/POLSPEVP3', 'POLYOMAVIRUS/POLSPEST', 'POLYOMAVIRUS/POLSPELT', 'RHINOVIRUS/RHISPECG', 'DENGE/DENSPECG', 'INFLUENZA/INSUBFNA', 'INFLUENZA/INFSUBHA', 'INFLUENZA/INFSUBMP', 'EBOLA/EBOSPECG']

iterations_number = 100
#dataset_path = "/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/"
dataset_path = sys.argv[1]


for i, dataset in enumerate(datasets):
    print(i, "\n\nEVALUATING DATASET: ", dataset, "...")
    print('===========================================================')
    print('===========================================================')    
    trainingData = fe.generateLabeledData(dataset_path + dataset + "/data.fa", dataset_path  + dataset + "/class.csv")

    f = open(dataset_path + dataset + "/results.metrics", "w")
    f.write("dataset;acc;precision;recall;fscore;nfeatures;k-mers;times\n")
    
    for j in range(iterations_number):  
        print('\niteration ....', j) 
        print('==================================')       
        shuffle(trainingData)
        train = trainingData[:int(len(trainingData)*0.8)]
        test = trainingData[int(len(trainingData)*0.8):len(trainingData)]
        #print(len(trainingData), len(train), len(test))
        
        best_k_mers, best_k_length, times = fe.getBestKmersAndFeatures('', train )
        print(dataset + "," + str(best_k_length) + "," + str(len(best_k_mers)) + "," + str(best_k_mers))

        print('training...')
        model = fe.train_model(train, best_k_mers)

        acc, precision, recall, fscore = fe.evaluation(model, test, best_k_mers)

        f.write(dataset + ";")
        f.write(str(acc) + ";")
        f.write(str(precision) + ";")
        f.write(str(recall) + ";")
        f.write(str(fscore) + ";")
        f.write(str(best_k_length) + ";")
        f.write(str(len(best_k_mers)) + ";")
        f.write(str(best_k_mers) + ";")
        f.write(str(times) + "\n")

    f.close()