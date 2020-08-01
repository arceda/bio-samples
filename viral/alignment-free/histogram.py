
# este script obtiene el skewness , kurtosis, energy de una sequencia de DNA , segun el paper: 
# DNA sequence similarity analysis using image texture analysis based on first-order statistics

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

current_dir = os.path.dirname(os.path.abspath(__file__))

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

def base_intensity(bases, G):
    if bases == "AA": intensity = G[0]
    if bases == "AG": intensity = G[1]
    if bases == "AC": intensity = G[2]
    if bases == "AT": intensity = G[3] 
    if bases == "GA": intensity = G[4] 
    if bases == "GG": intensity = G[5] 
    if bases == "GC": intensity = G[6] 
    if bases == "GT": intensity = G[7] 
    if bases == "CA": intensity = G[8]
    if bases == "CG": intensity = G[9] 
    if bases == "CC": intensity = G[10] 
    if bases == "CT": intensity = G[11]
    if bases == "TA": intensity = G[12]
    if bases == "TG": intensity = G[13] 
    if bases == "TC": intensity = G[14] 
    if bases == "TT": intensity = G[15] 
    return intensity

##############################################################################
# read fasta file 3
data = []
fa_file = "sample_genomes/KP317497.fna"
#fa_file = "sample_genomes/NR_117152.fna"
sequences = SeqIO.parse(fa_file, "fasta")
for record in sequences:
    data.append(record.seq.upper())
seq = data[0]

##############################################################################
############################      bases to numbers      ######################
##############################################################################
G = np.array([1,  17,  34,  51,  68,  85, 102, 119, 136, 153, 170, 187, 204, 221, 238, 255])

feature_vector = []
for i in range(0, len(seq)-1, 1):
    intensity = base_intensity(seq[i]+seq[i+1], G)    
    feature_vector.append( intensity )
feature_vector = np.array(feature_vector)
##############################################################################
##############################################################################


##############################################################################
############################    scale feature vector    ######################
##############################################################################
# 
'''
feature_vector = np.array(feature_vector)
c = feature_vector.min()
d = feature_vector.max()
feature_vector = feature_vector.astype(float) 
factor = (255 - 1)/(d-c)
feature_vector = ((feature_vector - c)*factor) + 1
'''
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
bins = np.append( [0], G+1 )

hist, bins = np.histogram(feature_vector, bins=bins)
hist = np.array(hist)

#df = pd.DataFrame({'x':G, 'y':hist})
#df.plot('x', 'y', kind='scatter')
#sns.distplot(hist, bins=bins, kde=False, rug=True);
sns.barplot(x=G, y=hist)
plt.savefig( current_dir + "/results/histogram.png" , dpi = 200, bbox_inches='tight')
##############################################################################
##############################################################################


##############################################################################
############################         features          #######################
##############################################################################
i_indices = G # used to avoid cycles, it represent the "i" in sum

p_hist = hist/hist.sum()  # equation (3)  p(i) = h(i)/NM

# equation (4)
mu = p_hist.dot( G.transpose() ) # equation (4)
print("mu", mu)
print(hist.mean(), p_hist.mean(), feature_vector.mean())
#mu = 0
#for index, intensity in enumerate(G):
#    mu += intensity*p_hist[index]
#print(mu)

# equation (5)
variance = (p_hist*(G - mu)**2).sum()
print("variance", variance, " variance numpy: ", np.var(feature_vector))
#variance_squared = 0
#for i variance range(G):
#   variance += (i - mu)*(i - mu)*p_hist[i]
#print(variance)

# equation (6)
standar_deviation = math.pow(variance, 1/2)
skewness = (1/math.pow(standar_deviation, 3)) * (p_hist*(G - mu)**3).sum()
print("skewness", skewness, " skew scipy:", skew(feature_vector)) #0,1027
#skewness = 0
#for i in range(G):
#    skewness += (i - mu)*(i - mu)*(i - mu)*p_hist[i]
#skewness = skewness/math.pow(standar_deviation, 3)
#print(skewness)

# equation (7)
my_kurtosis = (1/math.pow(standar_deviation, 4)) * ( (p_hist*(G - mu)**4) - 3 ).sum()
#my_kurtosis = (1/math.pow(standar_deviation, 4)) * ( (p_hist - 3)*(G - mu)**4 ).sum()
print("kurtosis", my_kurtosis, " kurtosis scipy:", kurtosis(feature_vector)) #-0,9176
#kurtosis = 0
#for i in range(G):
#    kurtosis += (i - mu)*(i - mu)*(i - mu)*(i - mu)*p_hist[i] - 3
#kurtosis = kurtosis/math.pow(standar_deviation, 4)
#print(kurtosis)

# equation (8)
energy = (p_hist*p_hist).sum()
print("energy", energy) #0,0739

# equation (9)
entropy = 0
for index, intensity in enumerate(G):
    if p_hist[index] != 0:
        entropy -= p_hist[index]*math.log( abs(p_hist[index]), 2 )
print("entropy", entropy) #3,9087

##############################################################################
##############################################################################









