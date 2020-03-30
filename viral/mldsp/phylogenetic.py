import numpy as np
from Bio import SeqIO
from Bio import Phylo

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

from skbio import DistanceMatrix
from skbio.tree import nj
from ete3 import PhyloTree, TreeStyle

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.cm as cm

from random import randint
from train import cmdscale

def readFasta(database):
    path = path_database + '/' + database

    number_of_clases = 0
    cluster_names = []
    points_per_cluster = []
    sequences = []
    str_all = ""
    
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
                str_all += ">" +   record.id + "\n"   +  str(record.seq.upper()) + "\n"       

    sequences_mat = np.array(sequences)

    return sequences_mat, number_of_clases, cluster_names, points_per_cluster, str_all

path_database = '/home/vicente/projects/BIOINFORMATICS/MLDSP/DataBase'
database_name = 'Dengue'
#path_database = sys.argv[1]
#database_name = sys.argv[2]

# readding database
print("readding model...")
current_dir = os.path.dirname(os.path.abspath(__file__))
file_name = current_dir + '/models/' + database_name
clf = joblib.load(file_name + ".sav")

print("readding magnitud spectrum...")
lg = np.genfromtxt(file_name + "_magnitud_spectrum.csv", delimiter=',', dtype=float)
print(lg.shape)
median_len = lg.shape[1]

sequences_mat, number_of_clases, cluster_names, points_per_cluster, str_all = readFasta(database_name)

# PEarson correlation coeficient
pearsoncorr = np.corrcoef(np.matrix(lg))
dist_mat = (1 - pearsoncorr)/2
rows, cols = dist_mat.shape

# check simetry
# some elements in diagonal are close to zero
# dist_mat = np.around(dist_mat, decimals=5)

for i in range(dist_mat.shape[0]):
    #if dist_mat[i][i] != 0.0:
    #    print(i, " ", dist_mat[i][i])

    dist_mat[i][i] = 0.0


asym = dist_mat - dist_mat.T
asym_ = np.where(asym != 0.0, 1, 0)
print(" no simetric distances:", np.sum(asym_))

dist_mat = dist_mat+ dist_mat.T - np.diag(np.diag(dist_mat))

asym = dist_mat - dist_mat.T
asym_ = np.where(asym != 0.0, 1, 0)
print(" no simetric distances:", np.sum(asym_))

#for i in range(rows):
#    for j in range(cols):





#########################################################################################################
# Classical multidimentional scaling
Y, evals = cmdscale(dist_mat)
print("dist_mat.shape:", dist_mat.shape)
print("Y.shape:", Y.shape)

ax = plt.axes(projection='3d')

# Data for a three-dimensional line
#zline = np.linspace(0, 15, 1000)
#xline = np.sin(zline)
#yline = np.cos(zline)
#ax.plot3D(xline, yline, zline, 'gray')

# Data for three-dimensional scattered points
zdata = Y[:,0]
xdata = Y[:,1]
ydata = Y[:,2]

#import itertools
#colors = itertools.cycle(["r", "b", "g"])
colors = cm.rainbow(np.linspace(0, 1, number_of_clases))


print(points_per_cluster)
print(xdata.shape)
tmp = 0
for i in range(number_of_clases):        
    ini = tmp
    end = ini + points_per_cluster[i]
    #print(ini, end)
    #ax.scatter3D(xdata[ini:end], ydata[ini:end], zdata[ini:end], alpha=0.6, c=next(colors))
    ax.scatter3D(xdata[ini:end], ydata[ini:end], zdata[ini:end], alpha=0.6, c=colors[i], label=cluster_names[i])
    tmp = end 
    
ax.legend()   

plt.show()





#########################################################################################################
# Phylogenetic tree


dm = DistanceMatrix(dist_mat, sequences_mat[:, 0])
tree = nj(dm)
newick_str = nj(dm, result_constructor=str)
t = PhyloTree(newick_str)
#t.link_to_alignment(alignment=str_all, alg_format="fasta")
t.show()