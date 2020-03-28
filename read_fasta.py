import numpy as np
from Bio import SeqIO
import re

fa_file = "/home/vicente/projects/BIOINFORMATICS/datasets/VIRAL/POLYOMAVIRUS/VP1.fas"

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
    