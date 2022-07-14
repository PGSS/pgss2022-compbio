#{} "defline" -> "species, group"
# Fragments_(whatever the enzyme is called) ie: Fragments_MseI -> {list of fragments}
#seq_length -> "int"

import csv
list = []
with open('./Sequence_Analyses.csv', 'r') as seqAFile:
    
    seqA = csv.reader(seqAFile, delimiter='|')
    counter = 0
    database = []
    
    next(seqA)
    for row in seqA:

        data = {}
        
        data["defline"] = row[0]
        data["seq_length"] = int(row[1])
        data["fragments_MseI"] = [int(x) for x in row[2].split(',')]
        data["fragments_Hpy188I"] = [int(x) for x in row[3].split(',')]

        database.append(data)
        
        next(seqA)
        