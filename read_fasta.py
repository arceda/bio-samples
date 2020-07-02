import numpy as np
from Bio import SeqIO
import re

fa_file = "P21333.fasta"

data = []
data_full = []

# Open the sequences file}
sequences = SeqIO.parse(fa_file, "fasta")


for record in sequences:
    data.append([record.id, record.seq.upper(), None])
    data_full.append([record])


temp = next(SeqIO.parse(fa_file, "fasta"))
print("%s %i" % (temp.id, len(temp)))

record_dict = SeqIO.index(fa_file, "fasta")
print(len(record_dict))
#print(record_dict.get_raw("AJ006022").decode())
#record_dict.close()
    

from Bio import SeqIO

sequences = SeqIO.parse("P21333.fasta", "fasta")
for record in sequences:
    data1 = str(record.seq.upper()) # the fasta file just have one sequence 

sequences = SeqIO.parse("Q8BTM8.fasta", "fasta")
for record in sequences:
    data2 = str(record.seq.upper()) # the fasta file just have one sequence  

print(data1)
print(data2)