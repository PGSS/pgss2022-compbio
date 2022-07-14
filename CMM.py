sample_exp =  { "seq_length" : 560,   "fragments_MseI" : [481,200],  "fragments_Hpy188I" : [607,428,264]}

sample_db = {"defline" : "S002165652 Lactobacillus delbrueckii subsp. null; E9D2-1	Lineage=Root;rootrank;Bacteria;domain;Firmicutes;phylum;Bacilli;class;Lactobacillales;order;Lactobacillaceae;family;Lactobacillus;genus",
 "seq_length": 844,
 "fragments_MseI": [824, 20],
 "fragments_Hpy188I": [337, 36, 177, 38, 256]}

def gel_match(gel_exp, gel_db):
    #compares overall length to database
    length_is_match = abs(gel_exp["seq_length"]-gel_db["seq_length"]) > 50
    mse_is_match = frags_match(gel_exp["fragments_MseI"], gel_db("fragments_Mse1"))
    hpy_is_match = frags_match(gel_exp["fragments_Hpy188I"], gel_db["fragments_Hpy188I"])
    return length_is_match and mse_is_match and hpy_is_match






def frags_match(frag_exp, frag_db):
    #compares MseI fragment digest to database
    #exclude close to full length
    #exclude <100 bp from db
    #count remaining fragments
    if 1==1:
    #if not frag_match(gelexp["fragments_MSeI"], gel_db["fragments_MSeI"]): <---dont know if this works or not
        return False
    return True



print(gel_match(sample_exp, sample_db))