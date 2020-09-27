# para usar: sudo apt-get install ncbi-blast+
# sudo pip3 install pyblastbio

from pyblast import BioBlast
from pyblast.utils import make_linear, make_circular
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from Bio import SeqIO
import time

queries = [
  SeqRecord(Seq("ACGTGATTCGTCGTGTAGTTGAGTGTTACGTTGCATGTCGTACGTGTGTAGTGTCGTGTAGTGCTGATGCTACGTGATCG"))
]
subjects = [
  SeqRecord(Seq("TCGTGTAGTTGAGTGTTACGTTGCATGTCGTTACGTGATCG"), id="aa1")
  ]

# pyblast requires a 'topology' annotation on the SeqRecords.
# we can make records circular or linear using `make_linear` or `make_circular` methods

t0 = time.time()

subjects = make_linear(subjects)
queries = make_linear(queries)

blast = BioBlast(subjects, queries)
results = blast.blastn()

t1 = time.time() - t0
#print(t1)
#print(results)


fa_file = "P21333.fasta"
sequences = SeqIO.parse(fa_file, "fasta")
queries = []
subjects = []
for record in sequences:
    queries.append( record )


#sequences = SeqIO.parse("viral/viral_classification/sample_genomes/HIV.B.fasta", "fasta")
sequences = SeqIO.parse("/home/vicente/datasets/MLDSP/Primates/seq/H1.txt", "fasta")
for record in sequences:
    subjects.append( record )


total_time = 0
for i in range(20):
    t0 = time.time()

    subjects = make_linear(subjects)
    queries = make_linear(queries)

    blast = BioBlast(subjects, queries)
    results = blast.blastn()

    t1 = time.time() - t0
    total_time += t1

print(total_time)
