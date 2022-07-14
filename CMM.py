import csv
#SAMPLE DATA BELOW

#replace with user input (also set +/- gap for seq len)
sample_exp =  { "seq_length" : 560,   "fragments_MseI" : [481,200],  "fragments_Hpy188I" : [607,428,264]}
sample_exp["fragments_MseI"] = [x for x in sample_exp["fragments_MseI"] if x <= sample_exp["seq_length"]]
sample_exp["fragments_Hpy188I"] = [x for x in sample_exp["fragments_Hpy188I"] if x <= sample_exp["seq_length"]]

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

#LENGTH MATCHING
def db_seqlen_filter(gel_exp, gel_db, buffer=100):
    db_lenfilt = []
    for bacteria in gel_db:
        if abs(gel_exp["seq_length"]-bacteria["seq_length"]) < buffer:
            db_lenfilt.append(bacteria)
    return db_lenfilt

#FRAGMENT MATCHING
def db_fragnum_filter(gel_exp, gel_db, buffer=100):
    #If more bands in experimental than digital = not a match 
    #If more bands in the digital than the experimental, may be a match
    db_fragfilt = []
    for bacteria in gel_db:
        fragments_msei = [x for x in bacteria["fragments_MseI"] if x >= buffer]
        fragments_hpy188i = [x for x in bacteria["fragments_Hpy188I"] if x >= buffer]
        if len(fragments_msei) >= len(gel_exp['fragments_MseI']) and len(fragments_hpy188i) >= len(gel_exp['fragments_Hpy188I']):
            db_fragfilt.append(bacteria)
    return db_fragfilt


    
db_init = db_import('./Sequence_Analyses.csv')
print(len(db_init))
db_lenfilt = db_seqlen_filter(sample_exp, db_init)
print(len(db_lenfilt))
db_fragfilt = db_fragnum_filter(sample_exp, db_lenfilt)
print(len(db_fragfilt))


