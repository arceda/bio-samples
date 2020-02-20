
import csv
import numpy as np
from Bio import SeqIO
import re
import sys

import matplotlib.pyplot as plt
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_recall_fscore_support
from sklearn.feature_selection import RFE
from sklearn.decomposition import TruncatedSVD
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold

from sklearn.svm import SVC
from sklearn.naive_bayes import MultinomialNB
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier

import time

###################################################################################################
#################################        READ DATASET            ##################################
# out: matrix [ [seq.id, seq, virus_type], ... ]
# fa_file: file with sequences FASTA
# cls_file: csv with id sequence and virus type
def generateLabeledData(fa_file, cls_file):    
    data = []

    # parse csv
    with open(cls_file) as f:
        reader = dict(csv.reader(f))

    # parse the sequences, it used GENERATOR (estos no es util para archivos grandes, no llevan todo a memoria y leen poco  a poco en cada invocación similara a los iteradores)
    for record in SeqIO.parse(fa_file, "fasta"):
        if record.id in reader:
            # Generate table [Id, Sequences, Class]
            data.append([record.id, record.seq.upper(), reader[record.id]])

    return data

def generateData(fa_file):
    data = []
   
    for record in SeqIO.parse(fa_file, "fasta"):
        data.append([record.id, record.seq.upper(), None])
   
    return data
###################################################################################################
###################################################################################################
    

###################################################################################################
#################################        GENERATE K-MERS         ##################################
# out: vector with all unique substring (k-mers) from training data
# data: trainign dataset [ [seq.id, seq, virus_type], ... ]
# k: k size in k-mers
def generate_K_mers(data, k):
    
	# List of k-mer
    K_mer = []
    dict = {}
    
    # Initialization of the dictionary
    # we generate unique k-mers for all training data
    for d in data:
        for i in range(0, len(d[1]) - k + 1, 1): 
            dict[d[1][i:i + k]] = 0
    		
    	# Remove patterns not use
    # in FASTA format there is other letter that 
    for key in dict.keys():
        if bool(re.match('^[ACGT]+$', str(key))) == True: K_mer.append(str(key))
    return K_mer
###################################################################################################
###################################################################################################


###################################################################################################
#################################   OCURRENCE OF EACH K-MER      ##################################
# IT IS USED IN PREDICTIONS
def generateMatrice(data, K_mer, k):
	# Variables
    X = []

    # Generate K-mer dictionnary
    X_dict = {}
    for i, e in enumerate(K_mer):  
        X_dict[e] = 0

    #print('K_mer', K_mer)
    #print('len X_dict', len(X_dict))
	
	# Generates X (matrix attributes)
    for d in data:
        x = []
        x_dict =  X_dict.copy()

        # Count K-mer occurences (with overlaping)
        for i in range(0, len(d[1]) - k + 1, 1):
            try: x_dict[d[1][i:i + k]] = x_dict[d[1][i:i + k]] + 1; 
            except: pass

        # Get only occurences from dictionnary
        for value in x_dict:
            x.append(x_dict.get(value))
        X.append(x)

	# Return matrices X (matrix attributes)
    return X

# k_mrs: vector with all unique substring (k-mers) from training data
# data: trainign dataset [ [seq.id, seq, virus_type], ... ]
# k: k size in k-mers
#out: matrix of ocurrences de each k-mer(substring) per sequences (row:seq, col:k-mer), y: labels
def generateXYMatrice(data, K_mer, k):
	# Variables
	X = generateMatrice(data, K_mer, k)
	y = []

	# Generates y (matrix class)
	for i in data:
		y.append(i[2])
	
	# Return matrices X and y  (matrix attributes and matrix class)
	return X, y
###################################################################################################
###################################################################################################


###################################################################################################
#################################   MAX MIN NORMALIZATION        ##################################
def maxMinNormalization(X):
    X_max = max(max(X))
    if X_max > 1:		
        #print("Apply linearly scaling each attribute to the range [0, 1]")
        minMaxScaler = MinMaxScaler(feature_range=(0, 1), copy = False)
        X = minMaxScaler.fit_transform(X)
    #else: print("Scaling not required ")
    
    return X

###################################################################################################
###################################################################################################
  
    
###################################################################################################
#############################     RECURSIVE FEATURE ELIMINATION     ###############################
def recursiveFeatureElimination(X, y, k_mers, features_max):
    preliminary_rfe_step = 0.1 #10%
    #preliminary_rfe_step = 1
    clf = SVC(kernel = "linear", C = 1)   
    
    # original
    
    if len(X[0]) > features_max:
        #print("Preliminary - RFE...")	
        rfe = RFE(estimator = clf, n_features_to_select = features_max, step = preliminary_rfe_step)
        new_X = rfe.fit_transform(X, y)

        # Update list of k_mers
        for i, value in enumerate(rfe.support_):
            if value == False: k_mers[i] = None
        new_k_mers = list(filter(lambda a: a != None, k_mers))

        print("reduce features to ", len(new_k_mers))

        return new_X, new_k_mers 

        
  
    else:
        return X, k_mers    
    

def SVD(X, nfeatures = None):
    #print('performing SVD')
    
    _X = np.matrix(X)
    print('_X.shape', _X.shape)
    #print('_X', _X)
    
    if nfeatures == None:
        mayor_zero = np.count_nonzero(_X > 0.0) # count the > 0 entries
        #print('mayor_zero', mayor_zero)
        mayor_zero = mayor_zero/len(X) # mean 
        nfeatures = int(mayor_zero*0.1) # take the 10%
    
    print('nfeaturesn SVD', nfeatures)
    
    svd =  TruncatedSVD(algorithm='randomized', n_components = nfeatures)
    #svd =  TruncatedSVD(n_components = 54)
    X_transf = svd.fit_transform(_X)
    

    print('X_transf.shape', X_transf.shape)
    #print('X_transf', X_transf)
    
    return X_transf, nfeatures
        
    
###################################################################################################
################################################################################################### 
    


###################################################################################################
#############################        EVALUATE FEATURE SIZES         ###############################
# X: ocurrence matrix of each k-mer
def evaluateFeatureSizes(X, y, k_mers, range_features, features_max, n_splits):
    clf = SVC(kernel = "linear", C = 1)  
    
    scores = []
    supports = []

    for n in range_features:
        print("\rRFE :", round(n / features_max * 100, 0), "%", end='')
        f_measure = 0
        k_mers_rfe = [] # here we store the k-mers no eliminated
        rfe = RFE(estimator = clf, n_features_to_select = n, step = 1)
        X_rfe = rfe.fit_transform(X, y)

        # Update list of k_mers
        for i, j in enumerate(rfe.support_):
            if j == True: k_mers_rfe.append(k_mers[i])
    
        # Evaluation of attributes with F-measure and Cross-validation
        for train_index, test_index in StratifiedKFold(n_splits = n_splits, shuffle=False, random_state=None).split(X_rfe, y):
            X_train, X_test = list(X_rfe[test_index]), list(X_rfe[train_index])
            y_train, y_test = list(y[test_index]), list(y[train_index])

            # Prediction

            #start_time = time.clock()
            clf.fit(X_train, y_train)
            #print(time.clock() - start_time, "seconds")
            y_pred = clf.predict(X_test)
        
            # Calcul metric scores
            f_measure = f_measure + f1_score(y_test, y_pred, average ="weighted")

        # Calcul mean F-measure
        mean_f_measure = f_measure / n_splits

        # Save scores
        scores.append(mean_f_measure)
        supports.append(k_mers_rfe)
        
    return scores, supports
    
    
###################################################################################################
################################################################################################### 
    


###################################################################################################
#############################        GET OPTIMAL SOLUTION           ###############################
def getOptimalSolution(scores_list, supports_list, k_mers_range, features_range, T):
    best_score = 0
	# Optimal score in relation with treshold
    optimal_score = 0
    
    # Identify best solution
    # here we indetified the best set of k-mers (one element in supports_list) based on the max score (scores_list) 
    for i, s in enumerate(scores_list):
        if max(s) > best_score:
            best_score = max(s)
            index = s.index(max(s))
            best_k_length = k_mers_range[i] #here we store the value of k in k-mer
            best_k_mers = supports_list[i][index]
        elif max(s) == best_score:
            if s.index(max(s)) < index:
                best_score = max(s)
                index = s.index(max(s))
                best_k_length = k_mers_range[i]
                best_k_mers = supports_list[i][index]
        else: pass

	# Identify optimal solution
    # here, we considered all the solutions (scores_list)  that at least have the T% of the best score 
    # and considered the one that have  the fewest features
    for i, l in enumerate(scores_list):
        for j, s in enumerate(l):
            if s >=  best_score * T and j <= index: 
                optimal_score = s
                index = j
                best_k_length = k_mers_range[i]
                best_k_mers = supports_list[i][index]
                print("\nChange optimal solution")
			
    if optimal_score == 0: optimal_score = best_score

	# Save plot results
    fig = plt.figure(figsize=(12, 10))
    for i, s in enumerate(scores_list):
        label = str(k_mers_range[i]) + "-mers"
        plt.plot(features_range, s, label= label)
    plt.ylabel("F-measure")
    plt.xlabel("Number of features")
    plt.axvline(index + 1, linestyle=':', color='r')
    title = "F-measure : " + str(optimal_score) + " K-mer size : " + str(best_k_length) + " Number of features : " + str(index + 1)
    plt.title(title)
    plt.legend()
   
    
    fname = str('fig_file')
    plt.savefig(fname)

    return best_k_mers, best_k_length
 
###################################################################################################
###################################################################################################  
# it return the best k-mers and nfeatures based on the paper Toward an Alignment-Free Method for Feature Extraction
# and Accurate Classification of Viral Sequences
def getBestKmersAndFeatures(path, trainingData=None):
    
    
    features_max = 100
    features_min = 1
    n_splits = 5
    k_min = 1
    k_max = 30  
    T = 0.99

    range_k_mers = range(k_min, k_max + 1, 1)
    range_features = range(features_min, features_max + 1, 1)
    '''
    features_max = 20
    features_min = 5
    n_splits = 5
    k_min = 5
    k_max = 15  
    T = 0.99

    range_k_mers = range(k_min, k_max + 1, 5)
    range_features = range(features_min, features_max + 1, 5)
    '''

    scores_list = []
    supports_list = []
    
    if trainingData == None:
        trainingData = generateLabeledData(path + "/data.fa", path +  "/class.csv")
    #trainingData = generateLabeledData("../castor_krfe/Data/HIVGRPCG/data.fa", "../castor_krfe/Data/HIVGRPCG/class.csv")
    #data         = generateData("../castor_krfe/Data/HIVGRPCG/data.fa")
    
    start_time_t = time.clock()

    for k in range_k_mers: 
        print("\n\n Evaluating with k-mer:", k)
        
        start_time = time.clock()
        k_mers      = generate_K_mers(trainingData, k) #list of substring of size k: (if k = 2; k_mers= [AT, CG, AC, ...])    
        t_1 = time.clock() - start_time; print('generate_K_mers took', t_1, "seconds", "feature size ", len(k_mers))

        start_time = time.clock()
        X, y        = generateXYMatrice(trainingData, k_mers, k) # OCURERNCE MATRIX
        t_2 = time.clock() - start_time; print('generateXYMatrice took', t_2, "seconds")

        start_time = time.clock()
        X           = maxMinNormalization(X)
        t_3 = time.clock() - start_time; print('maxMinNormalization took', t_3, "seconds")

        # THIS TKE TOO LONG TIME!!!!!, WE START WITH 1K FEATURES
        start_time = time.clock()
        #print("feature size ", len(k_mers))
        X, k_mers   = recursiveFeatureElimination(X, y, k_mers, features_max)
        t_4 = time.clock() - start_time; print('recursiveFeatureElimination took', t_4, "seconds")
                           
    
        labelEncodel = LabelEncoder()
        y = labelEncodel.fit_transform(y)
    
        # score = f1-measure, support = list of k-mers
        start_time = time.clock()
        scores, supports =  evaluateFeatureSizes(X, y, k_mers, range_features, features_max, n_splits)
        t_5 = time.clock() - start_time; print('evaluateFeatureSizes took', t_5, "seconds")
    
        scores_list.append(scores)
        supports_list.append(supports) 


    start_time = time.clock()
    best_k_mers, best_k_length = getOptimalSolution(scores_list, supports_list, range_k_mers, range_features, T)
    t_6 = time.clock() - start_time; print('getOptimalSolution took', t_6, "seconds")

    t_0 = time.clock() - start_time_t

    return best_k_mers, best_k_length, [t_0, t_1, t_2, t_3, t_4, t_5, t_6]
   
###################################################################################################
#############################               TRAINING                ###############################

def train_model(training_data, best_k_mers):

    # Generate  matrices
    best_k_length = len(best_k_mers[0])
    X_train, y_train = generateXYMatrice(training_data, best_k_mers, best_k_length)

    # Implement and fit classifier
    clf = SVC(kernel = "linear", C = 1) 
    clf.fit(X_train, y_train)
    
    return clf

def evaluation(clf, testing_data, best_k_mers):

    # Generate matrices
    best_k_length   = len(best_k_mers[0])
    X_test, y_test  = generateXYMatrice(testing_data, best_k_mers, best_k_length)


    # Realize prediction
    y_pred = clf.predict(X_test)
        
    # Calcul metric scores
    metrics = precision_recall_fscore_support(y_test, y_pred, average='weighted')
    acc = accuracy_score(y_test, y_pred)

    print('metrics: acc, precision, recall, fscore ', acc, metrics)
    return acc, metrics[0], metrics[1], metrics[2]

###################################################################################################
################################################################################################### 






if __name__ == "__main__" :
    # CMD
    # python3 viral/viral_classification/feature_extractor.py '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/PAPILLOMA/HPVSPECG'


    #training_data = generateLabeledData("../Data/HIVGRPCG/data.fa", "../Data/HIVGRPCG/class.csv")
    path = sys.argv[1] # folder path if fasta and csv
    #path = '/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/PAPILLOMA/HPVSPECG'
    best_k_mers, best_k_length, times = getBestKmersAndFeatures(path)
    print("Identified k =", best_k_length)
    print("Identified k-mers  =", best_k_mers)

  