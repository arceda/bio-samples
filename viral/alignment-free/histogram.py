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
import math

def contrast_stretching(img):
    a = 0
    b = 255
    c = img.min()
    d = img.max()
    img = img.astype(float) 
    
    factor = (b-a)/(d-c)
    new_img = ((img - c)*factor) + a   

    for i, row in enumerate(new_img):
        for j, x in enumerate(row):
            if x < 0:
                new_img[i][j] = 0
            if x > 255:
                new_img[i][j] = 255

    return new_img

def base_intensity(bases):
    if bases == "AA": intensity = 0 
    if bases == "AG": intensity = 1 
    if bases == "AC": intensity = 2
    if bases == "AT": intensity = 3 
    if bases == "GA": intensity = 4 
    if bases == "GG": intensity = 5 
    if bases == "GC": intensity = 6 
    if bases == "GT": intensity = 7 
    if bases == "CA": intensity = 8 
    if bases == "CG": intensity = 9 
    if bases == "CC": intensity = 10 
    if bases == "CT": intensity = 11
    if bases == "TA": intensity = 12
    if bases == "TG": intensity = 13 
    if bases == "TC": intensity = 14 
    if bases == "TT": intensity = 15 
    return intensity

##############################################################################
# read fasta file 3
data = []
fa_file = "sample_genomes/KP317497.fna"
sequences = SeqIO.parse(fa_file, "fasta")
for record in sequences:
    data.append(record.seq.upper())
seq = data[0]

##############################################################################
############################      bases to numbers      ######################
##############################################################################

feature_vector = []
for i in range(0, len(seq)-1, 1):
    intensity = base_intensity(seq[i]+seq[i+1])    
    feature_vector.append( intensity )

##############################################################################
##############################################################################


##############################################################################
############################    scale feature vector    ######################
##############################################################################
# 
feature_vector = np.array(feature_vector)
c = feature_vector.min()
d = feature_vector.max()
feature_vector = feature_vector.astype(float) 
factor = (255 - 1)/(d-c)
feature_vector = ((feature_vector - c)*factor) + 1

##############################################################################
##############################################################################



##############################################################################
# vector to image
width = 70
height = math.ceil( len(seq)/70 )
feature_vector_img = np.append(feature_vector, np.zeros( ( width*height - len(feature_vector) , 1) )) # completas con ceros
#print(feature_vector.shape)
img = feature_vector_img.reshape((height, width))
#print(img.shape) # rows, cols

img = img.astype('int8')

cv2.imwrite("gene.jpg", img)
#cv2.imshow("win", img)
#cv2.waitKey()

##############################################################################
############################        histogram          #######################
##############################################################################
G = 256
hist, bins = np.histogram(feature_vector, bins=range(0,G+1))
hist = np.array(hist)
##############################################################################
##############################################################################


##############################################################################
############################         features          #######################
##############################################################################
i_indices = np.array(range(0, G)) # used to avoid cycles, it represent the "i" in sum

p_hist = hist/feature_vector.shape[0]  # equation (3)  p(i) = h(i)/NM

# equation (4)
mu = p_hist.dot( i_indices.transpose() ) # equation (4)
print(mu)
#mu = 0
#for i in range(G):
#    mu += i*p_hist[i]
#print(mu)

# equation (5)
variance_squared = ((i_indices - mu)*(i_indices - mu)*p_hist).sum()
print(variance_squared)
#variance_squared = 0
#for i in range(G):
#   variance_squared += (i - mu)*(i - mu)*p_hist[i]
#print(variance_squared)

# equation (6)
standar_deviation = math.pow(variance_squared, 1/2)
skewness = (1/math.pow(standar_deviation, 3)) * ((i_indices - mu)*(i_indices - mu)*(i_indices - mu)*p_hist).sum()
print(skewness)
#skewness = 0
#for i in range(G):
#    skewness += (i - mu)*(i - mu)*(i - mu)*p_hist[i]
#skewness = skewness/math.pow(standar_deviation, 3)
#print(skewness)

# equation (7)
kurtosis = (1/math.pow(standar_deviation, 4)) * ((i_indices - mu)*(i_indices - mu)*(i_indices - mu)*(i_indices - mu)*p_hist - 3).sum()
print(kurtosis)
#kurtosis = 0
#for i in range(G):
#    kurtosis += (i - mu)*(i - mu)*(i - mu)*(i - mu)*p_hist[i] - 3
#kurtosis = kurtosis/math.pow(standar_deviation, 4)
#print(kurtosis)

# equation (8)
energy = (p_hist*p_hist).sum()
print(energy)

# equation (9)
entropy = 0
for i in range(G):
    if p_hist[i] != 0:
        entropy += p_hist[i]*math.log( abs(p_hist[i]), 2 )
print(entropy)

##############################################################################
##############################################################################









