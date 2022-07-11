def is_match(gel_exp, gel_db):
    #compares overall length to expected length
    if abs(gel_exp["seq_length"]-gel_db["seq_length"]) > 50:
        return False
    return True