import csv
#SAMPLE DATA BELOW

#replace with user input (also set +/- gap for seq len)
sample_exp =  { "seq_length" : 560,   "fragments_MseI" : [481,200],  "fragments_Hpy188I" : [607,428,264]}

sample_db = {"defline" : "S002165652 Lactobacillus delbrueckii subsp. null; E9D2-1	Lineage=Root;rootrank;Bacteria;domain;Firmicutes;phylum;Bacilli;class;Lactobacillales;order;Lactobacillaceae;family;Lactobacillus;genus",
"seq_length": 844,
"fragments_MseI": [824, 20],
"fragments_Hpy188I": [337, 36, 177, 38, 256]}




#GEL MATCHING

def db_cut(gel_exp, gel_db, buffer=50):
    db = []
    for bacteria in gel_db:
        if abs(gel_exp["seq_length"]-bacteria["seq_length"]) < buffer:
            db.append(bacteria)
    return db


def gel_match(gel_exp, gel_db, buffer=50):
    #compares overall length to database
    length_is_match = abs(gel_exp["seq_length"]-gel_db["seq_length"]) > buffer
    mse_is_match = frags_match(gel_exp["fragments_MseI"], gel_db["fragments_MseI"])
    hpy_is_match = frags_match(gel_exp["fragments_Hpy188I"], gel_db["fragments_Hpy188I"])
    return length_is_match and mse_is_match and hpy_is_match



#FRAGMENT MATCHING


# def frags_match(frag_exp, frag_db):
#     #compares MseI fragment digest to database
#     #exclude close to full length
#     #exclude <100 bp from db
#     #count remaining fragments
#     if 1==1:
#     #if not frags_match(gelexp["fragments_MSeI"], gel_db["fragments_MSeI"]): 
#         return False
#     return True

#     if gel_exp["fragments_MSEI"]



#FORMATS SEQUENCE LIST

def db_import(db_path):
    with open(db_path, 'r') as seqFile:
        
        #creates list of definitions for all sequences
        seq = csv.reader(seqFile, delimiter='|')
        database = []
        

        #iterates through sequence list
        next(seq)
        for row in seq:

            data = {}
            
            data["defline"] = row[0]
            data["seq_length"] = int(row[1])
            data["fragments_MseI"] = [int(x) for x in row[2].split(',')]
            data["fragments_Hpy188I"] = [int(x) for x in row[3].split(',')]

            #adds the definitions to a list
            database.append(data)
            
            next(seq)
    return database

db_init = db_import('./Sequence_Analyses.csv')
db = db_cut(sample_exp, db_init)
print(len(db))
#databs = db_import(db_init)