
import csv
import numpy as np
from Bio import SeqIO
import re
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
def generateMatrice(data, K_mer, k):
	# Variables
    X = []

    # Generate K-mer dictionnary
    X_dict = {}
    for i, e in enumerate(K_mer):  
        X_dict[e] = 0
	
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
        print("Apply linearly scaling each attribute to the range [0, 1]")
        minMaxScaler = MinMaxScaler(feature_range=(0, 1), copy = False)
        X = minMaxScaler.fit_transform(X)
    else: print("Scaling not required ")
    
    return X

###################################################################################################
###################################################################################################
  
    
###################################################################################################
#############################     RECURSIVE FEATURE ELIMINATION     ###############################
def recursiveFeatureElimination(X, y, k_mers, features_max):
    preliminary_rfe_step = 0.1
    clf = SVC(kernel = "linear", C = 1)   
    
    if len(X[0]) > features_max:
        print("Preliminary - RFE...")	
        rfe = RFE(estimator = clf, n_features_to_select = features_max, step = preliminary_rfe_step)
        new_X = rfe.fit_transform(X, y)

        # Update list of k_mers
        for i, value in enumerate(rfe.support_):
            if value == False: k_mers[i] = None
        new_k_mers = list(filter(lambda a: a != None, k_mers))

    return new_X, new_k_mers    
    
    
###################################################################################################
################################################################################################### 
    


###################################################################################################
#############################        EVALUATE FEATURE SIZES         ###############################
def evaluateFeatureSizes(X, y, k_mers, features_min, features_max):
    clf = SVC(kernel = "linear", C = 1)   

    for n in range(features_min, features_max + 1):
        print("\rRFE :", round(n / features_max * 100, 0), "%", end='')
        f_measure = 0
        k_mers_rfe = []
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
            clf.fit(X_train, y_train)
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
    for i, s in enumerate(scores_list):
        if max(s) > best_score:
            best_score = max(s)
            index = s.index(max(s))
            best_k_length = k_mers_range[i]
            best_k_mers = supports_list[i][index]
        elif max(s) == best_score:
            if s.index(max(s)) < index:
                best_score = max(s)
                index = s.index(max(s))
                best_k_length = k_mers_range[i]
                best_k_mers = supports_list[i][index]
        else: pass

	# Identify optimal solution
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
 
###################################################################################################
###################################################################################################  
    

if __name__ == "__main__" :
    #training_data = generateLabeledData("../Data/HIVGRPCG/data.fa", "../Data/HIVGRPCG/class.csv")
    
    features_max = 50
    features_min = 10
    n_splits = 5
    k_min = 3
    k_max = 4  
    T = 0.1

    scores_list = []
    supports_list = []
    
    trainingData = generateLabeledData("../castor_krfe/Data/HIVGRPCG/data.fa", "../castor_krfe/Data/HIVGRPCG/class.csv")
    data         = generateData("../castor_krfe/Data/HIVGRPCG/data.fa")
    
    for k in range(k_min, k_max + 1): 
        
        k_mers      = generate_K_mers(trainingData, k)    
        X, y        = generateXYMatrice(trainingData, k_mers, k) # OCURERNCE MATRIX
        X           = maxMinNormalization(X)
        X, k_mers   = recursiveFeatureElimination(X, y, k_mers, features_max)
                
    
        scores = []
        supports = []
    
        labelEncodel = LabelEncoder()
        y = labelEncodel.fit_transform(y)
    
        scores, supports =  evaluateFeatureSizes(X, y, k_mers, features_min, features_max)
    
        scores_list.append(scores)
        supports_list.append(supports)


    getOptimalSolution(scores_list, supports_list, range(k_min, k_max + 1), range(features_min, features_max + 1), T)
   

    '''
    fa_file = "../castor_krfe/Data/HIVGRPCG/data.fa"
    cls_file = "../castor_krfe/Data/HIVGRPCG/class.csv"
    
    data = []

    # Open the class file
    with open(cls_file) as f:
        reader = dict(csv.reader(f))

    # Open the sequences file}
    sequences = SeqIO.parse(fa_file, "fasta")


    for record in sequences:
        
        
        
        if record.id in reader:
            # Generate table [Id, Sequences, Class]
            data.append([record.id, record.seq.upper(), reader[record.id]])


    temp = next(SeqIO.parse("../castor_krfe/Data/HIVGRPCG/data.fa", "fasta"))
    print("%s %i" % (temp.id, len(temp)))
    
    record_dict = SeqIO.index("../castor_krfe/Data/HIVGRPCG/data.fa", "fasta")
    print(len(record_dict))
    print(record_dict.get_raw("AJ006022").decode())
    record_dict.close()
    '''