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
# the bases to numbers
feature_vector = []
for i in range(0, len(seq)-1, 1):
    intensity = base_intensity(seq[i]+seq[i+1])    
    feature_vector.append( intensity )

##############################################################################
# scale feature vector
feature_vector = np.array(feature_vector)
c = feature_vector.min()
d = feature_vector.max()
feature_vector = feature_vector.astype(float) 
factor = (255)/(d-c)
feature_vector = ((feature_vector - c)*factor)   
#print(feature_vector)

##############################################################################
# vector to image
width = 70
height = math.ceil( len(seq)/70 )
feature_vector = np.append(feature_vector, np.zeros( ( width*height - len(feature_vector) , 1) )) # completas con ceros
#print(feature_vector.shape)
img = feature_vector.reshape((height, width))
#print(img.shape) # rows, cols

img = img.astype('int8')

cv2.imwrite("gene.jpg", img)
#cv2.imshow("win", img)
#cv2.waitKey()

##############################################################################
# histogram
hist, bins = np.histogram(feature_vector, bins=range(0,257))
hist = np.array(hist)

print(hist.shape)
mu = hist.dot( np.array(range(0,256)).transpose())
print(mu)
