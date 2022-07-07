from fasta_reader import read_fasta
from Bio.Seq import Seq
from Bio.Restriction import *
rb = RestrictionBatch(['MseI', 'Hpy188I'])

for item in read_fasta("C:/Users/AWESD/Downloads/lactic_acid_only.fa.gz"): #change this to wherever you have the file
    print(item)
    seq = Seq(item.sequence)
    a = Analysis(rb, seq)
    a.print_that()
    break #remove this to print everything in the database