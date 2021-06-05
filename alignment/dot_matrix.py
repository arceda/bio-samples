import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import pylab

def delta(x,y):
    return 0 if x == y else 1

def M(seq1,seq2,i,j,k):
    return sum(delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))

def makeMatrix(seq1,seq2,k):
    n = len(seq1)
    m = len(seq2)
    return [[M(seq1,seq2,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]

def plotMatrix(M,t, seq1, seq2, nonblank = chr(0x25A0), blank = ' '):
    print(' |' + seq2)
    print('-'*(2 + len(seq2)))
    for label,row in zip(seq1,M):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line)

def dotplot(seq1,seq2,k = 1,t = 1):
    M = makeMatrix(seq1,seq2,k)
    plotMatrix(M, t, seq1,seq2) #experiment with character choice

# read sequences
sequences = SeqIO.parse("P21333.fasta", "fasta")
for record in sequences:
    data1 = str(record.seq.upper()) # the fasta file just have one sequence 

sequences = SeqIO.parse("Q8BTM8.fasta", "fasta")
for record in sequences:
    data2 = str(record.seq.upper()) # the fasta file just have one sequence  

#print(data1)
#print(data2)

data1 = "ACCTGATACAGTGGCT"
data2 = "ACCTGATAGATACAGTGGCT"

# simple dot matrix in console
#dotplot(data1,data2)

#dotplot=plt.imshow(np.array(makeMatrix(data1,data2,1)))
#xt=plt.xticks(np.arange(len(list(data1))),list(data1))
#yt=plt.yticks(np.arange(len(list(data1))),list(data1))
#plt.show()

# from biopython

window = 3
seq_one = data1
seq_two = data2
data = [
    [
        (seq_one[i : i + window] != seq_two[j : j + window])
        for j in range(len(seq_one) - window)
    ]
    for i in range(len(seq_two) - window)
]

pylab.gray()
pylab.imshow(data)
pylab.xlabel("%s (length %i bp)" % ("seq1", len(seq_one)))
pylab.ylabel("%s (length %i bp)" % ("seq2", len(seq_two)))
#pylab.title("Dot plot using window size %i\n(allowing no mis-matches)" % window)
pylab.title("Dot plot using window size %i" % window)
pylab.show()

#from biopython with scatterplot
dict_one = {}
dict_two = {}
for (seq, section_dict) in [
    (str(data1), dict_one),
    (str(data2), dict_two),
]:
    for i in range(len(seq) - window):
        section = seq[i : i + window]
        try:
            section_dict[section].append(i)
        except KeyError:
            section_dict[section] = [i]
# Now find any sub-sequences found in both sequences
# (Python 2.3 would require slightly different code here)
matches = set(dict_one).intersection(dict_two)
print("%i unique matches" % len(matches))

# Create lists of x and y co-ordinates for scatter plot
x = []
y = []
for section in matches:
    for i in dict_one[section]:
        for j in dict_two[section]:
            x.append(i)
            y.append(j)


pylab.cla()  # clear any prior graph
pylab.gray()
pylab.scatter(x, y)
pylab.xlim(0, len(data1) - window)
pylab.ylim(0, len(data2) - window)
pylab.xlabel("%s (length %i bp)" % ("seq1", len(data1)))
pylab.ylabel("%s (length %i bp)" % ("seq2", len(data2)))
#pylab.title("Dot plot using window size %i\n(allowing no mis-matches)" % window)
pylab.title("Dot plot using window size %i" % window)
pylab.show()