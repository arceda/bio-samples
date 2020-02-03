# Imports
from sklearn.svm import SVC
from sklearn.naive_bayes import MultinomialNB
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier

# SVM
def svm(): 
    return SVC(kernel = "linear", C = 1)

# Naive Bayes
def multinomialNB(): 
    return MultinomialNB(alpha = 1.0, fit_prior = True, class_prior = None)

# KNN
def kNeighborsClassifier(): 
    return KNeighborsClassifier(n_neighbors = 5, algorithm = "auto", metric = "minkowski", n_jobs = 4)

# Random Forest
def randomForestClassifier(): 
    return RandomForestClassifier(n_estimators = 10, criterion = "gini", max_depth = None, bootstrap = True, n_jobs = 4, random_state = 0)

# Multi-layer Perceptron
def mLPClassifier(): 
    return MLPClassifier(hidden_layer_sizes=(100, ), activation = "relu", solver = "adam", alpha = 0.0001, learning_rate="constant", learning_rate_init = 0.001, max_iter = 200,  random_state = 0, tol = 0.0001)

