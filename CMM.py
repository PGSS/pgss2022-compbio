import csv
#SAMPLE DATA BELOW

#replace with user input (also set +/- gap for seq len)
sample_exp =  { "seq_length" : 560,   "fragments_MseI" : [481,200],  "fragments_Hpy188I" : [607,428,264]}
sample_exp["fragments_MseI"] = [x for x in sample_exp["fragments_MseI"] if x <= sample_exp["seq_length"]]
sample_exp["fragments_Hpy188I"] = [x for x in sample_exp["fragments_Hpy188I"] if x <= sample_exp["seq_length"]]

#LENGTH MATCHING
def db_len_filter(gel_exp, gel_db, buffer=100):
    db = []
    for bacteria in gel_db:
        if abs(gel_exp["seq_length"]-bacteria["seq_length"]) < buffer:
            db.append(bacteria)
    return db





#FRAGMENT MATCHING
def frags_sort(frag_exp, frag_db):
    #compares MseI fragment digest to database
    #exclude close to full length
    #exclude <100 bp from db
    #count remaining fragments
    if 1==1:
    #if not frags_match(gelexp["fragments_MSeI"], gel_db["fragments_MSeI"]): 
        return False
    return True

    if gel_exp["fragments_MSEI"]

def 



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
print(len(db_init))
db = db_len_filter(sample_exp, db_init)
print(len(db))

