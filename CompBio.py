def is_length(gel_exp, gel_db):
    #compares overall length to database
    if abs(gel_exp["seq_length"]-gel_db["seq_length"]) > 50:
        return False
    return True
def is_frag_msei(frag_exp, frag_db):
    #compares MseI fragment digest to database
    if 1 == 1:
        return False
    return True
def is_frag_hpy188i(frag_exp, frag_db):
    #compares Hpy188I fragment digest to database
    if 1 == 1:
        return False
    return True
