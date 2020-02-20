# this script compare the performace of:
# Toward an Alignment-Free Method for Feature Extraction [1]
# An open-source k-mer based machine learning tool for fast and accurate [2]

import matplotlib.pyplot as plt
from sklearn.metrics import f1_score
from sklearn.feature_selection import RFE
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_recall_fscore_support

from sklearn.svm import SVC
from sklearn.naive_bayes import MultinomialNB
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold
import time
import feature_extractor as fe
from random import shuffle
import sys
import numpy as np

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


###################################################################################################
#############################               PREOCESS FEATURES       ###############################
# this funtion return de acc, fcore, .... of the models proposed by
# Toward an Alignment-Free Method for Feature Extraction [1]
# An open-source k-mer based machine learning tool for fast and accurate [2]
# they differs only in the feature reduction tecnhique.
def getFeatures(train, test, k, nfeatures, clf):
    start_time = time.clock()
    k_mers      = fe.generate_K_mers(train, k) #list of substring of size k: (if k = 2; k_mers= [AT, CG, AC, ...])    
    t_1 = time.clock() - start_time; print('generate_K_mers took', t_1, "seconds", "feature size ", len(k_mers))

    _k_mers = k_mers.copy()

    start_time = time.clock()
    X, y        = fe.generateXYMatrice(train, k_mers, k) # OCURRENCE MATRIX
    t_2 = time.clock() - start_time; print('generateXYMatrice took', t_2, "seconds")

    start_time = time.clock()
    X           = fe.maxMinNormalization(X)
    t_3 = time.clock() - start_time; print('maxMinNormalization took', t_3, "seconds")

    # dimentionality reduction with recursiveFeatureElimination [1]
    start_time = time.clock()
    X_1, k_mers_1   = fe.recursiveFeatureElimination(X, y, k_mers, nfeatures)
    t_4 = time.clock() - start_time; print('recursiveFeatureElimination took', t_4, "seconds")

    # dimentionality reduction with recursiveFeatureElimination [2]
    start_time = time.clock()
    X_2, nfeatures_svd  = fe.SVD(X)
    t_5 = time.clock() - start_time; print('SVD took', t_5, "seconds")
                        

    #labelEncodel = LabelEncoder()
    #y = labelEncodel.fit_transform(y)


    # train with recursiveFeatureElimination [1]
    ############################################
    print('training...[1]') #in this alg no perfomrm max min normalization in training
    model = fe.train_model(train, k_mers_1)
    acc, precision, recall, fscore = fe.evaluation(model, test, k_mers_1)

    # train with recursiveFeatureElimination [2]
    ############################################
    print('training...[2]')   
    '''
    clf = SVC(kernel = "linear", C = 1) 
    clf.fit(X_2, y)  
    X_test, y_test          = fe.generateXYMatrice(test, _k_mers, k)
    #print('k_mers', _k_mers)  
    #print('nfeatures_svd', nfeatures_svd)   
    X_test                  = fe.maxMinNormalization(X_test)
    X_test                  = fe.SVD(X_test, nfeatures_svd)          
    
    y_pred = clf.predict(X_test)
    
    metrics = precision_recall_fscore_support(y_test, y_pred, average='weighted')
    acc = accuracy_score(y_test, y_pred)
    print(acc, metrics)
'''
    

###################################################################################################
###################################################################################################

def split_data(data, indexs):
    k_fold = []
    for i in indexs:
        k_fold.append(data[i])
    
    return k_fold

k = 7
nfeatures = 7

datasets = ['POLYOMAVIRUS/POLSPEVP1', 'POLYOMAVIRUS/POLSPEVP2', 'POLYOMAVIRUS/POLSPEVP3', 
            'POLYOMAVIRUS/POLSPEST', 'POLYOMAVIRUS/POLSPELT'] 

dataset_path = "/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/"
#dataset_path = sys.argv[1]

clf = SVC(kernel = "linear", C = 1)  


for i, dataset in enumerate(datasets):

    print(i, "\n\nEVALUATING DATASET: ", dataset, "...")
    print('===========================================================')
    print('===========================================================') 

    for k in range(5,11):
        print("\n\n Evaluating with k-mer:", k)

        data = fe.generateLabeledData(dataset_path + dataset + "/data.fa", dataset_path  + dataset + "/class.csv")         
        
        kf = KFold(n_splits=5, shuffle=True, random_state=1)
        for train_index, test_index in kf.split(np.zeros(len(data))):
            print("\n\n Evaluating k-fold")

            data_train = split_data(data, train_index)
            data_test = split_data(data, test_index)                       

            getFeatures(data_train, data_test, k, nfeatures, clf)
            
            

        
    break

