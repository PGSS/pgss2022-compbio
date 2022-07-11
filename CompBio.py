sample_exp =  { "seq_length" : 560,   "fragments_MseI" : [481,200],  "fragments_Hpy188I" : [607, 428,264]}

sample_db = {"defline" : "S002165652 Lactobacillus delbrueckii subsp. null; E9D2-1	Lineage=Root;rootrank;Bacteria;domain;Firmicutes;phylum;Bacilli;class;Lactobacillales;order;Lactobacillaceae;family;Lactobacillus;genus",
 "seq_length": 844,
 "fragments_MseI": [824, 20],
 "fragments_Hpy188I": [337, 36, 177, 38, 256]}

def gel_match(gel_exp, gel_db):
    #compares overall length to database
    if abs(gel_exp["seq_length"]-gel_db["seq_length"]) > 50:
        return False
    return True
def frags_match(frag_exp, frag_db):
    #compares MseI fragment digest to database
    #exclude close to full length
    #exclude <100 bp from db
    #count remaining fragments
    if 1 == 1:
        return False
    return True



print(gel_match(sample_exp, sample_db))