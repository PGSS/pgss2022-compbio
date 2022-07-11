import csv
from typing import Sequence
from Bio.Seq import Seq
from Bio.Restriction import *
from fasta_reader import read_fasta

analyzed_list = []
rb = RestrictionBatch(['MseI', 'Hpy188I'])




def get_fragment_from(enzyme_cut, length):
  output = ""
  x = 0
  for i in enzyme_cut:
    output += str(i - x)
    output += ", "
    x = i
  output += str(length - x)
  return output


#with open('./Sequence_Analyses.csv', 'w') as csvfile:
    # creating a csv writer object 
    #csvwriter = csv.writer(csvfile)

h = open('./Sequence_Analyses.csv', 'w')
writer = csv.writer(h, delimiter='|')


processed = 0
a = None
for item in read_fasta("./lactic_acid_only.fa.gz"):
    bacteria_def = item.defline
    #print(bacteria_def)

    Sequence_Length = len(item.sequence)
    #print(f"Length of the sequence : {Sequence_Length}")

    analyzed_list.append(Analysis(rb.search(Seq(item.sequence)), Seq(item.sequence)))
    a = Analysis(rb.search(Seq(item.sequence)), Seq(item.sequence))
    #a.print_that()

    results = a.full()
    MseI_cut = results[MseI]
    MseI_lengths = get_fragment_from(MseI_cut, Sequence_Length)

    Hpy188I_cut = results[Hpy188I]
    Hpy188I_lengths = get_fragment_from(Hpy188I_cut, Sequence_Length)
    #print(MseI_lengths)
    #print(Hpy188I_lengths)

    writer.writerow([bacteria_def, Sequence_Length, MseI_lengths, Hpy188I_lengths,])


    processed += 1
    if processed % 100 == 0 :   
        print("Completed", processed)
    
    

h.close()
