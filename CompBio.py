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
