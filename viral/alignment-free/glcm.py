
# este es una replica del paper Use of image texture analysis to find DNA sequence similarities

from sklearn.model_selection import KFold 
from sklearn.model_selection import train_test_split
from matplotlib import pyplot as plt
from matplotlib import cm
import math
import os
import sys
import cv2
import numpy as np
import math
from scipy.stats import kurtosis, skew
from Bio import SeqIO
import pandas as pd
import seaborn as sns


import numpy as np
from skimage.feature import greycomatrix, greycoprops
from skimage import io, color, img_as_ubyte

from descriptor import get_features_glcm

##############################################################################
# read fasta file 3
data = []
#fa_file = "sample_genomes/KP317497.fna"
fa_file = "sample_genomes/L00016.fna"
sequences = SeqIO.parse(fa_file, "fasta")
for record in sequences:
    data.append(record.seq.upper())
seq = data[0]

##############################################################################
############################      bases to numbers      ######################
##############################################################################

feature_vector = []
for i, base in enumerate(seq):
    if base == 'A': val = 1
    if base == 'C': val = 2
    if base == 'G': val = 3
    if base == 'T': val = 4

    intensity = val + i
    feature_vector.append( intensity )
feature_vector = np.array([feature_vector], dtype=np.uint8)
levels = np.max(feature_vector) + 1
#print(feature_vector.shape)
#print(feature_vector)
#print(levels)
##############################################################################
##############################################################################

##############################################################################
############################             GLCM           ######################
##############################################################################

glcm_full = greycomatrix(feature_vector, [1], [0], levels=levels, symmetric=False, normed=True)
#glcm_full = greycomatrix(feature_vector, [1], [0, np.pi/4, np.pi/2, 3*np.pi/4], levels=levels, symmetric=False, normed=False)
glcm = glcm_full[:, :, 0, 0]
#print(glcm)
#print(np.sum(glcm))

# https://scikit-image.org/docs/stable/api/skimage.feature.html#skimage.feature.greycoprops
# contrast  ###############################################3
contrast = greycoprops(glcm_full, 'contrast')
my_contrast = 0
for i in range(glcm.shape[0]):
    for j in range(glcm.shape[1]):
        my_contrast += ((i-j)**2)*glcm[i][j]
print('contrast ', contrast, my_contrast)

# correlation ###############################################
correlation = greycoprops(glcm_full, 'correlation')
print('correlation ', correlation)

# energy ###############################################
energy = greycoprops(glcm_full, 'energy')
print('energy ', energy)

# homogeneity ###############################################
homogeneity = greycoprops(glcm_full, 'homogeneity')
print('homogeneity ', homogeneity)


# contrast  ###############################################3
entropy = -(glcm*np.ma.log(np.abs(glcm))).sum()

my_entropy = 0
for i in range(glcm.shape[0]):
    for j in range(glcm.shape[1]):
        if glcm[i][j] != 0:
            my_entropy += math.log( abs(glcm[i][j]))*glcm[i][j]

print('entropy ', entropy, my_entropy)


##############################################################################
############################             GLCM           ######################
##############################################################################

# GLCM example
image = np.array([[0, 0, 1, 1],
                   [0, 0, 1, 1],
                   [0, 2, 2, 2],
                   [2, 2, 3, 3]], dtype=np.uint8)

image = np.array([[0, 0, 1, 1, 2, 2, 2, 3]], dtype=np.uint8)

#compute GLCM with 4 different angles, levels are the max intensity and the number of rows cols
result = greycomatrix(image, [1], [0, np.pi/4, np.pi/2, 3*np.pi/4], levels=4)

#compute GLCM horizontal, levels are the max intensity and the number of rows cols
result = greycomatrix(image, [1], [0], levels=4)
#print(result[:, :, 0, 0])
result = greycomatrix(image, [1], [0, np.pi/4, np.pi/2, 3*np.pi/4], levels=4)
#print(result[:, :, 0, 0])


if __name__ == "__main__" :
    data = []
    #fa_file = "sample_genomes/KP317497.fna"
    fa_file = "sample_genomes/V00672.fna"
    sequences = SeqIO.parse(fa_file, "fasta")
    for record in sequences:
        data.append(record.seq.upper())
    seq = data[0]

    print(get_features_glcm(seq))

    current_dir = os.path.dirname(os.path.abspath(__file__))
    sequences = [   'L00016.fna',       'M22650.fna',           'M22651.fna',           'M22653.fna',       'M22654.fna', 
                'M22655.fna',       'M22656.fna',           'M22657.fna',           'V00658.fna',       'V00659.fna', 
                'V00672.fna',       'V00675.fna']

    names    = [    'Human',            'Macaca mulatta',       'Macaca fuscata',       'Macaca fascicularis',  'Macaca sylvanus', 
                'Saimiri sciureus', 'Tarsius syrichta',     'Lemur catta',          'Gorilla',              'Hylobates', 
                'Chimpanzee',       'Sumatran Orangutan']

    data_features_glcm = []



    for sequence_file in sequences:       

        data = []    
        fa_file = current_dir + "/sample_genomes/" + sequence_file
        seqs = SeqIO.parse(fa_file, "fasta")
        for record in seqs:
            #print(record.id)
            #print(record.description)
            data.append(record.seq.upper())      

        seq = data[0]           

        entropy, contrast, energy, correlation, homogeneity = get_features_glcm(seq)
        data_features_glcm.append( [entropy, contrast, energy, correlation, homogeneity] )

    data_features_glcm = np.array(data_features_glcm)

    
    x = data_features_glcm[:, 2]
    y = data_features_glcm[:, 0]


    plt.plot(x, y, 'o', color='black')

    for i, txt in enumerate(names):
        plt.annotate(txt, (x[i], y[i]), textcoords="offset points", xytext=(-15,-15), ha='left')


    plt.show()