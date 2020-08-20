# this script compute Chaos game representation of DNA
# extracted from https://towardsdatascience.com/chaos-game-representation-of-a-genetic-sequence-4681f1a67e14


import collections
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm
import pylab
import math
import os
import sys
import cv2
import numpy as np

from Bio import SeqIO
import re

#if not sys.warnoptions:
#    import warnings
#    warnings.simplefilter("ignore")
 
def count_kmers(sequence, k):
    d = collections.defaultdict(int)
    for i in range(len(sequence)-(k-1)):
        d[sequence[i:i+k]] +=1

    #print(d.keys)
    #for key in d.keys():

        #print(key)
        #if "N" in key:
        #    del d[key]
    return d
 
def probabilities(sequence, kmer_count, k):
    probabilities = collections.defaultdict(float)
    N = len(sequence)
    for key, value in kmer_count.items():
        probabilities[key] = float(value) / (N - k + 1)
    return probabilities
 
def chaos_game_representation(probabilities, k):
    array_size = int(math.sqrt(4**k))
    chaos = []
    for i in range(array_size):
        chaos.append([0]*array_size)
 
    maxx = array_size
    maxy = array_size
    posx = 1
    posy = 1
    for key, value in probabilities.items():
        for char in key:
            if char == "T":
                posx += maxx / 2                
            elif char == "C":
                posy += maxy / 2
            elif char == "G":
                posx += maxx / 2
                posy += maxy / 2
            maxx = maxx / 2
            maxy /= 2
        chaos[int(posy-1)][int(posx-1)] = value
        maxx = array_size
        maxy = array_size
        posx = 1
        posy = 1
 
    return chaos

def save_fcgr(id, sequence, k, path, name = ''):
    chaos = chaos_game_representation(probabilities(str(sequence), count_kmers(str(sequence), k), k), k)
    
    # show with pylab
    #pylab.title('Chaos game representation for ' + str(k) + '-mers')
    #pylab.imshow(chaos, interpolation='nearest', cmap=cm.gray_r)
    #pylab.show()
    #pylab.savefig(current_dir + '/cgr_k=' + str(k) + '.png', dpi = 300)

    #save with opencv
    img = np.matrix(chaos)
    min_value = img.min()
    max_value = img.max()
    img = ((img - min_value)/(max_value-min_value))*255
    #img = 255 - img  # the same al matplot
    #print(img.shape, img)    
    #cv2.imwrite(path + '/' + name + '_' + id + '_k=' + str(k) + '.jpg', img)
    cv2.imwrite(path + '/' + name + '_k=' + str(k) + '.jpg', img)
    
 
if __name__ == "__main__" :

    current_dir = os.path.dirname(os.path.abspath(__file__))
    path_sequence = sys.argv[1]
    k = int(sys.argv[2])
    #path_sequence = current_dir + "/sample-genomes/HIV.A1.fasta"
    #k = 8

    data = []
    sequences = SeqIO.parse(path_sequence, "fasta")
    for record in sequences:
        data.append([record.id, record.seq.upper()])

    #print(data)

    for seq in data:
        save_fcgr(seq[0], seq[1], k, current_dir)
    

