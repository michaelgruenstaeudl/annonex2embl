# Input: .nex, .csv

# 1. Parse data from NEXUS-file
from Bio.Nexus import Nexus
aln = Nexus.Nexus()
aln.read('/home/michael_science/Desktop/test.nex')
annos = aln.charsets
align = aln.matrix

# 2. Create SeqRecord
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#simple_seq = Seq("GATC")
for seqname in align.keys():
    seq = align[seqname]
    seqrec = SeqRecord(seq, id=seqname)
    seqrec.annotations = annos

    # 3. Convert SeqRecord to SeqFeature while degapping sequence
