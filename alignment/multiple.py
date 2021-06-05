from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

a = SeqRecord(Seq("AAAACGT", generic_dna), id="Alpha")
b = SeqRecord(Seq("AAA-CGT", generic_dna), id="Beta")
c = SeqRecord(Seq("AAAAGGT", generic_dna), id="Gamma")
align = MultipleSeqAlignment([a, b, c],
                             annotations={"tool": "demo"},
                             column_annotations={"stats": "CCCXCCC"})
print(align)
print(align.annotations)
print(align.column_annotations)