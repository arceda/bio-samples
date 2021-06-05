from Bio import pairwise2
from Bio.pairwise2 import format_alignment

"""
The match parameters are:

CODE  DESCRIPTION
x     No parameters. Identical characters have score of 1, otherwise 0.
m     A match score is the score of identical chars, otherwise mismatch
      score.
d     A dictionary returns the score of any pair of characters.
c     A callback function returns scores.
The gap penalty parameters are:

CODE  DESCRIPTION
x     No gap penalties.
s     Same open and extend gap penalties for both sequences.
d     The sequences have different open and extend gap penalties.
c     A callback function returns the gap penalties.
"""


alignments = pairwise2.align.globalxx("ACCGT", "ACG")
print(format_alignment(*alignments[0]))

alignments = pairwise2.align.localxx("ACCGT", "ACG")
print(format_alignment(*alignments[0]))
