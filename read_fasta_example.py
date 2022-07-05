from fasta_reader import read_fasta


for item in read_fasta("./lactic_acid_only.fa.gz"): #change this to wherever you have the file
    print(item.defline)
    break #remove this to print everything in the database
