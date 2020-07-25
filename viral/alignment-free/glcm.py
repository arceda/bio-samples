
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


##############################################################################
# read fasta file 3
data = []
#fa_file = "sample_genomes/KP317497.fna"
fa_file = "sample_genomes/NR_117152.fna"
sequences = SeqIO.parse(fa_file, "fasta")
for record in sequences:
    data.append(record.seq.upper())
seq = data[0]

##############################################################################
############################      bases to numbers      ######################
##############################################################################
G = np.array([1,  17,  34,  51,  68,  85, 102, 119, 136, 153, 170, 187, 204, 221, 238, 255])

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

