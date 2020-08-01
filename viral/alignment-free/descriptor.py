# Obtiene los atributos, segun los papers:
# [1] DNA sequence similarity analysis using image texture analysis based on first-order statistics
# [2] Use of image texture analysis to find DNA sequence similarities
# [3] A signal processing method for alignment-free metagenomic binning: multi-resolution genomic binary patterns
# 

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


def bases2number(seq, option):
    feature_vector = []
    if option == 0: #first order statistics [1]. each pair of bases are a number        
        # intensities
        G = np.array([1,  17,  34,  51,  68,  85, 102, 119, 136, 153, 170, 187, 204, 221, 238, 255])
        for i in range(0, len(seq)-1, 1):
            intensity = base_intensity(seq[i]+seq[i+1], G)    
            feature_vector.append( intensity )

    if option == 1: #used in glcm [2], each bases is a number + the i position
        for i, base in enumerate(seq):
            if base == 'A': val = 1
            if base == 'C': val = 2
            if base == 'G': val = 3
            if base == 'T': val = 4
            intensity = val + i
            feature_vector.append( intensity )

    if option == 2: #used in mlbp. Integer [3]
        for i, base in enumerate(seq):
            if base == 'A': val = 2
            if base == 'C': val = -1
            if base == 'G': val = 1
            if base == 'T': val = -2
            feature_vector.append( val )

    if option == 3: #used in mlbp. EIIP [3]
        for i, base in enumerate(seq):
            if base == 'A': val = 0.1260
            if base == 'C': val = 0.1340
            if base == 'G': val = 0.0806
            if base == 'T': val = 0.1335
            feature_vector.append( val )

    if option == 43: #used in mlbp. Atomic [3]
        for i, base in enumerate(seq):
            if base == 'A': val = 70
            if base == 'C': val = 58
            if base == 'G': val = 66
            if base == 'T': val = 78
            feature_vector.append( val )

    if option == 5: #used in mlbp. Real [3]
        for i, base in enumerate(seq):
            if base == 'A': val = -1.5
            if base == 'C': val = -0.5
            if base == 'G': val = 0.5
            if base == 'T': val = 1.5
            feature_vector.append( val )

    feature_vector = np.array(feature_vector)
    return feature_vector

def get_features(seq):
    # intensities
    G = np.array([1,  17,  34,  51,  68,  85, 102, 119, 136, 153, 170, 187, 204, 221, 238, 255])

    # feature vector
    feature_vector = bases2number(seq, 0)
    #feature_vector = []
    #for i in range(0, len(seq)-1, 1):
    #    intensity = base_intensity(seq[i]+seq[i+1], G)    
    #    feature_vector.append( intensity )
    #feature_vector = np.array(feature_vector)

    # histogram
    bins = np.append( [0], G+1 )
    hist, bins = np.histogram(feature_vector, bins=bins)
    hist = np.array(hist)
    p_hist = hist/hist.sum()

    mu = feature_vector.mean()
    variance = np.var(feature_vector)
    skewness = skew(feature_vector)
    my_kurtosis = kurtosis(feature_vector)
    energy = (p_hist*p_hist).sum()
    entropy = 0
    for index, intensity in enumerate(G):
        if p_hist[index] != 0:
            entropy -= p_hist[index]*math.log( abs(p_hist[index]), 2 )

    return skewness, my_kurtosis, energy, entropy


def get_features_glcm(seq):
    feature_vector = bases2number(seq, 1)
    #feature_vector = []
    #for i, base in enumerate(seq):
    ##    if base == 'A': val = 1
    #    if base == 'C': val = 2
    #    if base == 'G': val = 3
    #    if base == 'T': val = 4
    #
    #    intensity = val + i
    #    feature_vector.append( intensity )
    feature_vector = np.array([feature_vector], dtype=np.uint8)
    levels = np.max(feature_vector) + 1

    glcm_full = greycomatrix(feature_vector, [1], [0], levels=levels, symmetric=False, normed=True)
    glcm = glcm_full[:, :, 0, 0]

    entropy = -(glcm*np.ma.log(np.abs(glcm))).sum()
    contrast = greycoprops(glcm_full, 'contrast')    
    energy = greycoprops(glcm_full, 'energy')    
    correlation = greycoprops(glcm_full, 'correlation')        
    homogeneity = greycoprops(glcm_full, 'homogeneity')    
    
    return entropy, contrast[0][0], energy[0][0], correlation[0][0], homogeneity[0][0]

def sign(x):
    if x < 0 :
        return 0
    else:
        return 1

def get_features_lbp(seq):
    feature_vector = bases2number(seq, 2)
    #print(feature_vector.shape)
    p = 6
    lbp = []
    for t, x in enumerate(feature_vector):
        temp = 0
        for i in range( int(p/2) ):
            if t + i - p/2 >= 0 and t + i + 1 < feature_vector.shape[0]:
                temp = sign( feature_vector[t + i - int(p/2)] - feature_vector[t])*(2**i)  + sign( feature_vector[t + i + 1] - feature_vector[t])*(2**(i+int(p/2))) 

        lbp.append(temp)
    

    # histogram
    lbp_max = np.max(lbp)
    exponential = math.ceil(math.log(lbp_max, 2))
    bins = [2**i for i in range(exponential + 1)]
    #bins = np.append( [0], bins )

    #print(lbp_max)
    #print(exponential)
    #print(bins)

    hist, bins = np.histogram(lbp, bins=bins)
    #print(hist)
    #print(lbp)
    return np.array(hist)

def get_features_mlbp(seq):
    feature_vector = bases2number(seq, 2)
    #print(feature_vector.shape)
    p = 6
    result = []
    for z in range(2,p+1):
        lbp = []
        for t, x in enumerate(feature_vector):
            temp = 0
            p=z
            for i in range( int(p/2) ):
                if t + i - p/2 >= 0 and t + i + 1 < feature_vector.shape[0]:
                    temp = sign( feature_vector[t + i - int(p/2)] - feature_vector[t])*(2**i)  + sign( feature_vector[t + i + 1] - feature_vector[t])*(2**(i+int(p/2))) 

            lbp.append(temp)
        

        # histogram
        lbp_max = np.max(lbp)
        exponential = math.ceil(math.log(lbp_max, 2))
        bins = [2**i for i in range(exponential + 1)]
        hist, bins = np.histogram(lbp, bins=bins)
   
        result = np.hstack((result, hist))

    return np.array(result)

if __name__ == "__main__" :



    data = []
    #fa_file = "sample_genomes/KP317497.fna"
    fa_file = "sample_genomes/V00672.fna"
    sequences = SeqIO.parse(fa_file, "fasta")
    for record in sequences:
        data.append(record.seq.upper())
    seq = data[0]
    print(get_features_lbp(seq))
    print(get_features_mlbp(seq))
    