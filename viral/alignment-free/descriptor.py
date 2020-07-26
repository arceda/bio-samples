# Obtiene los atributos, segun el paper:
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

def get_features(seq):
    # intensities
    G = np.array([1,  17,  34,  51,  68,  85, 102, 119, 136, 153, 170, 187, 204, 221, 238, 255])

    # feature vector
    feature_vector = []
    for i in range(0, len(seq)-1, 1):
        intensity = base_intensity(seq[i]+seq[i+1], G)    
        feature_vector.append( intensity )
    feature_vector = np.array(feature_vector)

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

    glcm_full = greycomatrix(feature_vector, [1], [0], levels=levels, symmetric=False, normed=True)
    glcm = glcm_full[:, :, 0, 0]

    entropy = -(glcm*np.ma.log(np.abs(glcm))).sum()
    contrast = greycoprops(glcm_full, 'contrast')    
    energy = greycoprops(glcm_full, 'energy')    
    correlation = greycoprops(glcm_full, 'correlation')        
    homogeneity = greycoprops(glcm_full, 'homogeneity')    
    
    return entropy, contrast[0][0], energy[0][0], correlation[0][0], homogeneity[0][0]


if __name__ == "__main__" :
    data = []
    #fa_file = "sample_genomes/KP317497.fna"
    fa_file = "sample_genomes/V00672.fna"
    sequences = SeqIO.parse(fa_file, "fasta")
    for record in sequences:
        data.append(record.seq.upper())
    seq = data[0]

    print(get_features_glcm(seq))
    